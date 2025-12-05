use crate::{
    cli::ClusterAlgorithm,
    commands::shared,
    core::{MotifQuery, ReadRecord},
    features::extract::{
        self, AggregateRecord, ExtractionSink, FastqRecord, MotifSummaryRecord, PerReadRecord,
    },
    io::{writers::ExtractionWriter, BamReader, BamRecordDecoder, SequenceCache},
};
use anyhow::{anyhow, bail, Context, Result};
use bstr::ByteSlice;
use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use hdbscan::{Hdbscan, HdbscanHyperParams};
use noodles_bam::{self as bam, Record};
use noodles_bgzf as bgzf;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::{
    collections::{hash_map::Entry, HashMap, HashSet},
    f64::consts::PI,
    fs::{self, File},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Arc, Mutex,
    },
};

pub fn run(
    bam: PathBuf,
    motifs: Vec<String>,
    motif_file: Option<PathBuf>,
    sequence_fallback: Option<PathBuf>,
    sequence_index: Option<PathBuf>,
    output_dir: PathBuf,
    min_cluster_size: usize,
    min_samples: Option<usize>,
    emit_fastq: bool,
    threads: usize,
    cluster_algorithm: ClusterAlgorithm,
    bin_ids: Option<Vec<String>>,
    contigs: Vec<String>,
) -> Result<()> {
    if min_cluster_size == 0 {
        bail!("min-cluster-size must be greater than zero");
    }
    let motif_queries = shared::load_motif_queries(motifs, motif_file, bin_ids.as_deref())?;
    fs::create_dir_all(&output_dir)
        .with_context(|| format!("failed to create {}", output_dir.display()))?;

    let per_read_path = output_dir.join("per_read.tsv");
    let aggregate_path = output_dir.join("aggregate.tsv");
    let writer = ExtractionWriter::new(
        &motif_queries,
        per_read_path,
        aggregate_path,
        None::<PathBuf>,
        None::<PathBuf>,
    )?;
    eprintln!(
        "split_reads: extracting {} motifs from {}",
        motif_queries.len(),
        bam.display()
    );

    let mut reader =
        BamReader::open(&bam).with_context(|| format!("failed to open BAM {}", bam.display()))?;
    let sequence_cache = sequence_fallback
        .as_ref()
        .map(|path| {
            SequenceCache::from_path(path, sequence_index.as_deref())
                .with_context(|| format!("failed to load sequence fallback {}", path.display()))
        })
        .transpose()?
        .map(Arc::new);
    let decoder = BamRecordDecoder::new(reader.header(), sequence_cache);

    let contig_filter = if contigs.is_empty() {
        None
    } else {
        Some(contigs)
    };
    let mut sink_holder = Some(ClusterSink::new(writer, &motif_queries));
    let (sink, processed_reads) = if threads <= 1 {
        let mut extractor = crate::features::extract::Extractor::new(
            motif_queries.clone(),
            sink_holder.take().unwrap(),
        )?;
        let mut count = 0usize;
        reader.visit_records(contig_filter.as_ref().map(|c| c.as_slice()), |record| {
            let read = decoder.decode(record)?;
            extractor.process_read(&read)?;
            count += 1;
            if count % 1_000 == 0 {
                eprintln!("split_reads: processed {} reads...", count);
            }
            Ok(())
        })?;
        let sink = extractor.finish_with_sink()?;
        (sink, count)
    } else {
        let shared = SharedSink::new(sink_holder.take().unwrap());
        let extractor = ParallelExtractor::new(motif_queries.clone(), shared, threads)?;
        let count = extractor.run(
            &mut reader,
            &decoder,
            contig_filter.as_ref().map(|c| c.as_slice()),
        )?;
        let sink = extractor.into_sink()?;
        (sink, count)
    };
    decoder.report();
    eprintln!(
        "split_reads: finished extraction phase ({} reads); starting clustering",
        processed_reads
    );
    let (reads, raw_features, sequences, motif_names) = sink.into_data();
    let total_reads = reads.len();
    let (reads, raw_features, sequences) =
        drop_reads_without_features(reads, raw_features, sequences);
    let dropped_reads = total_reads.saturating_sub(reads.len());
    if dropped_reads > 0 {
        eprintln!(
            "split_reads: dropped {} reads with no motif calls before clustering",
            dropped_reads
        );
    }
    if reads.is_empty() {
        bail!("no reads matched motifs; nothing to cluster");
    }
    let (raw_features, motif_names, dropped_features) =
        retain_dense_features(raw_features, motif_names, 0.5)?;
    if dropped_features > 0 {
        eprintln!(
            "split_reads: dropped {} motifs with insufficient coverage",
            dropped_features
        );
    }
    if motif_names.is_empty() {
        bail!("all motifs were dropped after filtering for coverage");
    }
    eprintln!(
        "Imputing feature matrix for {} reads and {} motifs...",
        reads.len(),
        motif_names.len()
    );
    let imputed_features = impute_features(&raw_features);

    let effective_min_cluster_size = auto_min_cluster_size(reads.len());
    let effective_min_samples =
        min_samples.unwrap_or_else(|| auto_hdbscan_min_samples(reads.len()));
    eprintln!(
        "split_reads: clustering {} reads (min_cluster_size={}, min_samples={})",
        reads.len(),
        effective_min_cluster_size,
        effective_min_samples
    );
    let assignments = cluster_reads(
        &imputed_features,
        cluster_algorithm,
        effective_min_cluster_size,
        effective_min_samples,
    )?;
    write_clustering_results(
        &output_dir,
        &reads,
        &imputed_features,
        &motif_names,
        &assignments,
    )?;
    write_raw_clustering_results(
        &output_dir,
        &reads,
        &raw_features,
        &motif_names,
        &assignments,
    )?;
    write_cluster_summary(&output_dir, &motif_names, &assignments, &imputed_features)?;
    let clusters_dir = output_dir.join("clusters");
    if emit_fastq {
        write_cluster_fastqs(clusters_dir, &reads, &assignments, &sequences)?;
    } else {
        write_cluster_bams(&bam, clusters_dir, &reads, &assignments)?;
    }
    eprintln!(
        "split_reads: clustering complete; outputs written under {}",
        output_dir.display()
    );

    Ok(())
}

#[derive(Clone)]
struct SequenceInfo {
    sequence: Vec<u8>,
    qualities: Vec<u8>,
}

enum ExtractionEvent {
    PerRead(PerReadRecord<'static>),
    Aggregate(AggregateRecord<'static>),
    Summary(MotifSummaryRecord<'static>),
    Fastq(FastqRecord<'static>),
}

struct SharedSink<S> {
    sender: Sender<ExtractionEvent>,
    inner: Arc<Mutex<Option<S>>>,
    join_handle: Arc<Mutex<Option<std::thread::JoinHandle<()>>>>,
}

impl<S> Clone for SharedSink<S> {
    fn clone(&self) -> Self {
        Self {
            sender: self.sender.clone(),
            inner: Arc::clone(&self.inner),
            join_handle: Arc::clone(&self.join_handle),
        }
    }
}

impl<S> SharedSink<S>
where
    S: ExtractionSink + Send + 'static,
{
    fn new(inner: S) -> Self {
        let (tx, rx) = unbounded();
        let state = Arc::new(Mutex::new(Some(inner)));
        let state_clone = Arc::clone(&state);
        let join_handle = Arc::new(Mutex::new(None));
        let join_handle_clone = Arc::clone(&join_handle);
        let handle = std::thread::spawn(move || {
            if let Some(mut sink) = state_clone.lock().ok().and_then(|mut guard| guard.take()) {
                for event in rx.iter() {
                    match event {
                        ExtractionEvent::PerRead(record) => {
                            let _ = sink.record_per_read(record);
                        }
                        ExtractionEvent::Aggregate(record) => {
                            let _ = sink.record_aggregate(record);
                        }
                        ExtractionEvent::Summary(record) => {
                            let _ = sink.record_motif_summary(record);
                        }
                        ExtractionEvent::Fastq(record) => {
                            let _ = sink.record_fastq(record);
                        }
                    }
                }
                let _ = sink.finish();
                *state_clone.lock().unwrap() = Some(sink);
            }
        });
        *join_handle_clone.lock().unwrap() = Some(handle);
        Self {
            sender: tx,
            inner: state,
            join_handle,
        }
    }

    fn into_inner(self) -> Result<S> {
        drop(self.sender);
        if let Some(handle) = self.join_handle.lock().unwrap().take() {
            let _ = handle.join();
        }
        match Arc::try_unwrap(self.inner) {
            Ok(mutex) => mutex
                .into_inner()
                .map_err(|e| anyhow!("failed to unlock sink: {e}"))?
                .ok_or_else(|| anyhow!("sink not available")),
            Err(_) => Err(anyhow!("sink still has outstanding references at shutdown")),
        }
    }
}

impl<S: ExtractionSink + Send + 'static> ExtractionSink for SharedSink<S> {
    fn record_per_read(&mut self, record: PerReadRecord<'_>) -> Result<()> {
        self.sender
            .send(ExtractionEvent::PerRead(record.into_owned()))
            .map_err(|e| anyhow!(e.to_string()))
    }

    fn record_aggregate(&mut self, record: AggregateRecord<'_>) -> Result<()> {
        self.sender
            .send(ExtractionEvent::Aggregate(record.into_owned()))
            .map_err(|e| anyhow!(e.to_string()))
    }

    fn record_motif_summary(&mut self, record: MotifSummaryRecord<'_>) -> Result<()> {
        self.sender
            .send(ExtractionEvent::Summary(record.into_owned()))
            .map_err(|e| anyhow!(e.to_string()))
    }

    fn record_fastq(&mut self, record: FastqRecord<'_>) -> Result<()> {
        self.sender
            .send(ExtractionEvent::Fastq(record.into_owned()))
            .map_err(|e| anyhow!(e.to_string()))
    }

    fn finish(&mut self) -> Result<()> {
        Ok(())
    }
}

struct ParallelExtractor<S> {
    pool: ThreadPool,
    threads: usize,
    motifs: Vec<MotifQuery>,
    sink: SharedSink<S>,
}

impl<S> ParallelExtractor<S>
where
    S: ExtractionSink + Send + 'static,
{
    fn new(motifs: Vec<MotifQuery>, sink: SharedSink<S>, threads: usize) -> Result<Self> {
        let worker_count = threads.max(1);
        let pool = ThreadPoolBuilder::new()
            .num_threads(worker_count)
            .build()
            .context("failed to create thread pool for split_reads")?;
        Ok(Self {
            pool,
            threads: worker_count,
            motifs,
            sink,
        })
    }

    fn run(
        &self,
        reader: &mut BamReader,
        decoder: &BamRecordDecoder,
        contigs: Option<&[String]>,
    ) -> Result<usize> {
        let queue_depth = (self.threads * 4).max(1);
        let (tx, rx): (Sender<ReadRecord>, Receiver<ReadRecord>) = bounded(queue_depth);
        let processed = Arc::new(AtomicUsize::new(0));
        let aborted = Arc::new(AtomicBool::new(false));
        let error_slot: Arc<Mutex<Option<anyhow::Error>>> = Arc::new(Mutex::new(None));
        let (done_tx, done_rx) = bounded(self.threads);

        for _ in 0..self.threads {
            let rx = rx.clone();
            let motifs = self.motifs.clone();
            let mut sink = self.sink.clone();
            let processed = Arc::clone(&processed);
            let aborted = Arc::clone(&aborted);
            let error_slot = Arc::clone(&error_slot);
            let done_tx = done_tx.clone();
            self.pool.spawn(move || {
                let mut scratch = crate::features::extract::MotifScratch::new(motifs.len());
                while !aborted.load(Ordering::Relaxed) {
                    match rx.recv() {
                        Ok(read) => {
                            if let Err(err) = extract::process_read_with_motifs(
                                &read,
                                &motifs,
                                &mut sink,
                                &mut scratch,
                            ) {
                                aborted.store(true, Ordering::Relaxed);
                                Self::store_error(&error_slot, err);
                                break;
                            }
                            let total = processed.fetch_add(1, Ordering::Relaxed) + 1;
                            if total % 1_000 == 0 {
                                eprintln!("split_reads: processed {} reads...", total);
                            }
                        }
                        Err(_) => break,
                    }
                }
                let _ = done_tx.send(());
            });
        }

        drop(done_tx);
        drop(rx);

        if let Err(err) = reader.visit_records(contigs, |record| {
            if aborted.load(Ordering::Relaxed) {
                return Err(anyhow!("parallel extraction aborted"));
            }
            let read = decoder.decode(record)?;
            tx.send(read)
                .map_err(|send_err| anyhow!("failed to dispatch read: {send_err}"))?;
            Ok(())
        }) {
            aborted.store(true, Ordering::Relaxed);
            Self::store_error(&error_slot, err);
        }

        drop(tx);

        while done_rx.recv().is_ok() {}

        if let Some(err) = error_slot.lock().expect("error mutex poisoned").take() {
            return Err(err);
        }

        Ok(processed.load(Ordering::Relaxed))
    }

    fn into_sink(self) -> Result<S> {
        self.sink.into_inner()
    }

    fn store_error(slot: &Arc<Mutex<Option<anyhow::Error>>>, err: anyhow::Error) {
        if let Ok(mut guard) = slot.lock() {
            if guard.is_none() {
                *guard = Some(err);
            }
        }
    }
}

struct ClusterSink<W: ExtractionSink> {
    inner: W,
    motif_names: Vec<String>,
    motif_index: HashMap<String, usize>,
    read_features: HashMap<String, Vec<Option<f64>>>,
    sequences: HashMap<String, SequenceInfo>,
    hit_reads: HashSet<String>,
}

impl<W: ExtractionSink> ClusterSink<W> {
    fn new(inner: W, motifs: &[crate::core::MotifQuery]) -> Self {
        let motif_names: Vec<String> = motifs.iter().map(|m| m.raw().to_string()).collect();
        let motif_index = motif_names
            .iter()
            .enumerate()
            .map(|(idx, name)| (name.clone(), idx))
            .collect();
        Self {
            inner,
            motif_names,
            motif_index,
            read_features: HashMap::new(),
            sequences: HashMap::new(),
            hit_reads: HashSet::new(),
        }
    }

    fn into_data(
        self,
    ) -> (
        Vec<String>,
        Vec<Vec<Option<f64>>>,
        HashMap<String, SequenceInfo>,
        Vec<String>,
    ) {
        let mut read_ids: Vec<String> = self
            .read_features
            .keys()
            .filter(|id| self.hit_reads.contains(*id))
            .cloned()
            .collect();
        read_ids.sort();
        let features = read_ids
            .iter()
            .map(|id| self.read_features.get(id).cloned().unwrap_or_default())
            .collect();
        (read_ids, features, self.sequences, self.motif_names)
    }
}

impl<W: ExtractionSink> ExtractionSink for ClusterSink<W> {
    fn record_per_read(&mut self, record: PerReadRecord<'_>) -> Result<()> {
        self.hit_reads.insert(record.read_id.to_string());
        self.inner.record_per_read(record)
    }

    fn record_aggregate(&mut self, record: AggregateRecord<'_>) -> Result<()> {
        if let Some(idx) = self.motif_index.get(record.motif_name) {
            let entry = self
                .read_features
                .entry(record.read_id.to_string())
                .or_insert_with(|| vec![None; self.motif_names.len()]);
            entry[*idx] = record.mean_probability;
        }
        self.inner.record_aggregate(record)
    }

    fn record_motif_summary(&mut self, record: MotifSummaryRecord<'_>) -> Result<()> {
        self.inner.record_motif_summary(record)
    }

    fn record_fastq(&mut self, record: FastqRecord<'_>) -> Result<()> {
        self.sequences
            .entry(record.read_id.to_string())
            .or_insert_with(|| SequenceInfo {
                sequence: record.sequence.to_vec(),
                qualities: record.qualities.to_vec(),
            });
        self.inner.record_fastq(record)
    }

    fn finish(&mut self) -> Result<()> {
        self.inner.finish()
    }
}

fn cluster_reads(
    data: &[Vec<f64>],
    algorithm: ClusterAlgorithm,
    min_cluster_size: usize,
    min_samples: usize,
) -> Result<Vec<i32>> {
    match algorithm {
        ClusterAlgorithm::Hdbscan => cluster_with_hdbscan(data, min_cluster_size, min_samples),
        ClusterAlgorithm::Gmm => cluster_with_gmm(data, min_cluster_size),
        ClusterAlgorithm::Agglomerative => cluster_with_agglomerative(data, min_cluster_size),
    }
}

fn cluster_with_hdbscan(
    data: &[Vec<f64>],
    min_cluster_size: usize,
    min_samples: usize,
) -> Result<Vec<i32>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }
    let min_cluster_size = min_cluster_size.max(2);
    let min_samples = min_samples.max(1).min(min_cluster_size);
    if data.len() < min_cluster_size {
        return Ok(vec![-1; data.len()]);
    }
    let hyper_params = HdbscanHyperParams::builder()
        .min_cluster_size(min_cluster_size)
        .min_samples(min_samples)
        .build();
    let clustered_data = data.to_vec();
    let clusterer = Hdbscan::new(&clustered_data, hyper_params);
    clusterer
        .cluster()
        .map_err(|e| anyhow!("HDBSCAN clustering failed: {e}"))
}

fn cluster_with_gmm(data: &[Vec<f64>], min_cluster_size: usize) -> Result<Vec<i32>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }
    if data.len() == 1 {
        return Ok(vec![0]);
    }
    let max_components = determine_max_components(data.len(), min_cluster_size);
    if max_components < 2 {
        return Ok(vec![0; data.len()]);
    }
    let mut best_labels: Option<Vec<i32>> = None;
    let mut best_score = f64::NEG_INFINITY;
    let mut fallback: Option<Vec<i32>> = None;
    for components in 2..=max_components {
        let model = fit_gmm_diag(data, components, 100, 1e-3)?;
        let labels = assign_gmm_labels(&model, data);
        if fallback.is_none() {
            fallback = Some(labels.clone());
        }
        if let Some(score) = silhouette_score(data, &labels) {
            if score > best_score {
                best_score = score;
                best_labels = Some(labels);
            }
        }
    }
    if let Some(labels) = best_labels {
        Ok(labels)
    } else if let Some(labels) = fallback {
        Ok(labels)
    } else {
        Ok(vec![0; data.len()])
    }
}

fn cluster_with_agglomerative(data: &[Vec<f64>], min_cluster_size: usize) -> Result<Vec<i32>> {
    let n = data.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    if n == 1 {
        return Ok(vec![0]);
    }
    let mut target_clusters = determine_max_components(n, min_cluster_size);
    if target_clusters == 0 {
        target_clusters = 1;
    }
    target_clusters = target_clusters.min(n).max(1);
    if target_clusters <= 1 {
        return Ok(vec![0; n]);
    }
    let mst_edges = build_mst(data);
    Ok(assign_clusters_from_mst(n, &mst_edges, target_clusters))
}

fn determine_max_components(read_count: usize, min_cluster_size: usize) -> usize {
    if read_count <= 1 {
        return read_count;
    }
    let denom = min_cluster_size.max(1);
    let approx = (read_count / denom).max(2);
    approx.min(read_count).min(8)
}

fn silhouette_score(data: &[Vec<f64>], labels: &[i32]) -> Option<f64> {
    if data.len() < 2 {
        return None;
    }
    let mut clusters: HashMap<i32, Vec<usize>> = HashMap::new();
    for (idx, label) in labels.iter().enumerate() {
        clusters.entry(*label).or_default().push(idx);
    }
    if clusters.len() < 2 {
        return None;
    }
    let mut total = 0.0;
    for (idx, sample) in data.iter().enumerate() {
        let cluster_id = labels[idx];
        let members = clusters.get(&cluster_id)?;
        let a = if members.len() <= 1 {
            0.0
        } else {
            let mut sum = 0.0;
            for &other in members {
                if other == idx {
                    continue;
                }
                sum += euclidean_distance(sample, &data[other]);
            }
            sum / (members.len() as f64 - 1.0)
        };
        let mut b = f64::INFINITY;
        for (other_id, other_members) in clusters.iter() {
            if *other_id == cluster_id {
                continue;
            }
            let mut sum = 0.0;
            for &other in other_members {
                sum += euclidean_distance(sample, &data[other]);
            }
            let avg = sum / other_members.len() as f64;
            if avg < b {
                b = avg;
            }
        }
        let s = if a == 0.0 && b == 0.0 {
            0.0
        } else {
            (b - a) / b.max(a)
        };
        total += s;
    }
    Some(total / data.len() as f64)
}

fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (x, y) in a.iter().zip(b.iter()) {
        let diff = x - y;
        sum += diff * diff;
    }
    sum.sqrt()
}

#[derive(Clone)]
struct GmmComponent {
    weight: f64,
    mean: Vec<f64>,
    variance: Vec<f64>,
}

fn fit_gmm_diag(
    data: &[Vec<f64>],
    component_count: usize,
    max_iters: usize,
    tol: f64,
) -> Result<Vec<GmmComponent>> {
    if data.is_empty() {
        return Ok(Vec::new());
    }
    let dims = data[0].len();
    if dims == 0 {
        return Ok(vec![GmmComponent {
            weight: 1.0,
            mean: Vec::new(),
            variance: Vec::new(),
        }]);
    }
    let global_variance = compute_global_variance(data);
    let mut components = initialize_gmm_components(data, component_count, &global_variance);
    let mut prev_log_likelihood = f64::NEG_INFINITY;
    for _ in 0..max_iters {
        let (responsibilities, log_likelihood) = expectation_step(data, &components);
        let (updated_components, moved) =
            maximization_step(data, &responsibilities, &global_variance);
        components = updated_components;
        if (log_likelihood - prev_log_likelihood).abs() < tol && !moved {
            break;
        }
        prev_log_likelihood = log_likelihood;
    }
    Ok(components)
}

fn initialize_gmm_components(
    data: &[Vec<f64>],
    component_count: usize,
    global_variance: &[f64],
) -> Vec<GmmComponent> {
    let mut components = Vec::with_capacity(component_count);
    let n = data.len();
    if n == 0 {
        return components;
    }
    let step = (n / component_count).max(1);
    for comp_idx in 0..component_count {
        let index = (comp_idx * step) % n;
        let point = data[index].clone();
        components.push(GmmComponent {
            weight: 1.0 / component_count as f64,
            mean: point,
            variance: global_variance.to_vec(),
        });
    }
    components
}

fn expectation_step(data: &[Vec<f64>], components: &[GmmComponent]) -> (Vec<Vec<f64>>, f64) {
    let mut responsibilities = vec![vec![0.0; components.len()]; data.len()];
    let mut log_likelihood = 0.0;
    for (idx, sample) in data.iter().enumerate() {
        let mut log_probs = vec![0.0; components.len()];
        for (comp_idx, component) in components.iter().enumerate() {
            log_probs[comp_idx] = log_gaussian_diag(sample, component);
        }
        let log_sum = log_sum_exp(&log_probs);
        log_likelihood += log_sum;
        for comp_idx in 0..components.len() {
            responsibilities[idx][comp_idx] = (log_probs[comp_idx] - log_sum).exp();
        }
    }
    (responsibilities, log_likelihood)
}

fn maximization_step(
    data: &[Vec<f64>],
    responsibilities: &[Vec<f64>],
    global_variance: &[f64],
) -> (Vec<GmmComponent>, bool) {
    let n = data.len();
    let k = responsibilities.first().map(|r| r.len()).unwrap_or(0);
    let dims = data.first().map(|d| d.len()).unwrap_or(0);
    let mut nk = vec![0.0f64; k];
    let mut means = vec![vec![0.0f64; dims]; k];
    for (row_idx, sample) in data.iter().enumerate() {
        for comp_idx in 0..k {
            let r = responsibilities[row_idx][comp_idx];
            nk[comp_idx] += r;
            for dim in 0..dims {
                means[comp_idx][dim] += r * sample[dim];
            }
        }
    }

    let mut moved = false;
    for comp_idx in 0..k {
        if nk[comp_idx] > 1e-8 {
            for dim in 0..dims {
                means[comp_idx][dim] /= nk[comp_idx];
            }
        } else {
            moved = true;
            means[comp_idx] = data[comp_idx % data.len()].clone();
            nk[comp_idx] = 1e-6;
        }
    }

    let mut variances = vec![vec![0.0f64; dims]; k];
    for (row_idx, sample) in data.iter().enumerate() {
        for comp_idx in 0..k {
            let r = responsibilities[row_idx][comp_idx];
            for dim in 0..dims {
                let diff = sample[dim] - means[comp_idx][dim];
                variances[comp_idx][dim] += r * diff * diff;
            }
        }
    }

    let mut components = Vec::with_capacity(k);
    for comp_idx in 0..k {
        let mut component_variance = vec![0.0f64; dims];
        for dim in 0..dims {
            let mut value = variances[comp_idx][dim] / nk[comp_idx];
            if !value.is_finite() || value < 1e-6 {
                value = global_variance[dim];
            }
            component_variance[dim] = value;
        }
        let weight = (nk[comp_idx] / n as f64).max(1e-6);
        components.push(GmmComponent {
            weight,
            mean: means[comp_idx].clone(),
            variance: component_variance,
        });
    }
    (components, moved)
}

fn assign_gmm_labels(components: &[GmmComponent], data: &[Vec<f64>]) -> Vec<i32> {
    let mut labels = Vec::with_capacity(data.len());
    for sample in data {
        let mut best_idx = 0usize;
        let mut best_score = f64::NEG_INFINITY;
        for (comp_idx, component) in components.iter().enumerate() {
            let score = log_gaussian_diag(sample, component);
            if score > best_score {
                best_score = score;
                best_idx = comp_idx;
            }
        }
        labels.push(best_idx as i32);
    }
    labels
}

fn compute_global_variance(data: &[Vec<f64>]) -> Vec<f64> {
    let dims = data[0].len();
    let mut mean = vec![0.0f64; dims];
    for sample in data {
        for dim in 0..dims {
            mean[dim] += sample[dim];
        }
    }
    for dim in 0..dims {
        mean[dim] /= data.len() as f64;
    }
    let mut variance = vec![0.0f64; dims];
    for sample in data {
        for dim in 0..dims {
            let diff = sample[dim] - mean[dim];
            variance[dim] += diff * diff;
        }
    }
    for dim in 0..dims {
        variance[dim] = (variance[dim] / data.len() as f64).max(1e-6);
    }
    variance
}

fn log_gaussian_diag(sample: &[f64], component: &GmmComponent) -> f64 {
    let mut log_prob = component.weight.ln();
    let ln_norm = (2.0 * PI).ln();
    for dim in 0..sample.len() {
        let var = component.variance[dim].max(1e-6);
        let diff = sample[dim] - component.mean[dim];
        log_prob += -0.5 * ((diff * diff) / var + var.ln() + ln_norm);
    }
    log_prob
}

fn log_sum_exp(values: &[f64]) -> f64 {
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    let sum: f64 = values.iter().map(|value| (value - max).exp()).sum();
    max + sum.ln()
}

fn build_mst(data: &[Vec<f64>]) -> Vec<(usize, usize, f64)> {
    let n = data.len();
    let mut in_tree = vec![false; n];
    let mut min_dist = vec![f64::INFINITY; n];
    let mut parent = vec![usize::MAX; n];
    let mut edges = Vec::with_capacity(n.saturating_sub(1));
    in_tree[0] = true;
    for v in 1..n {
        min_dist[v] = euclidean_distance(&data[0], &data[v]);
        parent[v] = 0;
    }
    for _ in 1..n {
        let mut next = None;
        for v in 0..n {
            if !in_tree[v] && (next.is_none() || min_dist[v] < min_dist[next.unwrap()]) {
                next = Some(v);
            }
        }
        let Some(u) = next else {
            break;
        };
        in_tree[u] = true;
        edges.push((parent[u], u, min_dist[u]));
        for v in 0..n {
            if in_tree[v] || u == v {
                continue;
            }
            let dist = euclidean_distance(&data[u], &data[v]);
            if dist < min_dist[v] {
                min_dist[v] = dist;
                parent[v] = u;
            }
        }
    }
    edges
}

fn assign_clusters_from_mst(
    node_count: usize,
    edges: &[(usize, usize, f64)],
    cluster_count: usize,
) -> Vec<i32> {
    if node_count == 0 {
        return Vec::new();
    }
    if cluster_count <= 1 {
        return vec![0; node_count];
    }
    let mut sorted = edges.to_vec();
    sorted.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
    let skip = cluster_count.saturating_sub(1).min(sorted.len());
    let mut disjoint = DisjointSet::new(node_count);
    for (idx, (u, v, _)) in sorted.iter().enumerate() {
        if idx < skip {
            continue;
        }
        disjoint.union(*u, *v);
    }
    let mut cluster_map: HashMap<usize, i32> = HashMap::new();
    let mut next_id = 0i32;
    let mut labels = vec![0i32; node_count];
    for node in 0..node_count {
        let root = disjoint.find(node);
        let entry = cluster_map.entry(root).or_insert_with(|| {
            let id = next_id;
            next_id += 1;
            id
        });
        labels[node] = *entry;
    }
    labels
}

struct DisjointSet {
    parent: Vec<usize>,
    size: Vec<usize>,
}

impl DisjointSet {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            size: vec![1; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    fn union(&mut self, a: usize, b: usize) {
        let mut root_a = self.find(a);
        let mut root_b = self.find(b);
        if root_a == root_b {
            return;
        }
        if self.size[root_a] < self.size[root_b] {
            std::mem::swap(&mut root_a, &mut root_b);
        }
        self.parent[root_b] = root_a;
        self.size[root_a] += self.size[root_b];
    }
}
fn write_clustering_results(
    output_dir: &Path,
    read_ids: &[String],
    features: &[Vec<f64>],
    motifs: &[String],
    assignments: &[i32],
) -> Result<()> {
    let path = output_dir.join("read_clustering.tsv");
    let file =
        File::create(&path).with_context(|| format!("unable to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    write!(writer, "read_id\tcluster_id")?;
    for motif in motifs {
        write!(writer, "\t{}", motif)?;
    }
    writeln!(writer)?;

    for (idx, read_id) in read_ids.iter().enumerate() {
        write!(writer, "{}\t{}", read_id, assignments[idx])?;
        for value in &features[idx] {
            write!(writer, "\t{:.4}", value)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

fn write_raw_clustering_results(
    output_dir: &Path,
    read_ids: &[String],
    features: &[Vec<Option<f64>>],
    motifs: &[String],
    assignments: &[i32],
) -> Result<()> {
    let path = output_dir.join("read_clustering_raw.tsv");
    let file =
        File::create(&path).with_context(|| format!("unable to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    write!(writer, "read_id\tcluster_id")?;
    for motif in motifs {
        write!(writer, "\t{}", motif)?;
    }
    writeln!(writer)?;

    for (idx, read_id) in read_ids.iter().enumerate() {
        write!(writer, "{}\t{}", read_id, assignments[idx])?;
        for value in &features[idx] {
            match value {
                Some(v) => write!(writer, "\t{:.4}", v)?,
                None => write!(writer, "\tNA")?,
            }
        }
        writeln!(writer)?;
    }
    Ok(())
}

fn write_cluster_summary(
    output_dir: &Path,
    motifs: &[String],
    assignments: &[i32],
    features: &[Vec<f64>],
) -> Result<()> {
    let path = output_dir.join("cluster_summary.tsv");
    let file =
        File::create(&path).with_context(|| format!("unable to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    write!(writer, "cluster_id\tread_count")?;
    for motif in motifs {
        write!(writer, "\tcentroid_{}", motif)?;
    }
    writeln!(writer)?;

    let mut sums: HashMap<i32, (usize, Vec<f64>)> = HashMap::new();
    for (row_idx, cluster_id) in assignments.iter().enumerate() {
        let entry = sums
            .entry(*cluster_id)
            .or_insert_with(|| (0, vec![0.0f64; motifs.len()]));
        entry.0 += 1;
        for (col_idx, value) in features[row_idx].iter().enumerate() {
            entry.1[col_idx] += *value;
        }
    }

    let mut cluster_ids: Vec<i32> = sums.keys().copied().collect();
    cluster_ids.sort();
    for cluster_id in cluster_ids {
        let (total_reads, sum_vec) = sums.get(&cluster_id).unwrap();
        write!(writer, "{}\t{}", cluster_id, total_reads)?;
        for sum in sum_vec {
            write!(writer, "\t{:.4}", sum / *total_reads as f64)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn write_cluster_fastqs(
    clusters_dir: PathBuf,
    read_ids: &[String],
    assignments: &[i32],
    sequences: &HashMap<String, SequenceInfo>,
) -> Result<()> {
    fs::create_dir_all(&clusters_dir)
        .with_context(|| format!("failed to create {}", clusters_dir.display()))?;
    let mut writers: HashMap<i32, BufWriter<File>> = HashMap::new();

    let mut quality_buffer = Vec::new();
    for (idx, read_id) in read_ids.iter().enumerate() {
        let seq = match sequences.get(read_id) {
            Some(seq) => seq,
            None => continue,
        };
        let cluster_id = assignments[idx];
        if !writers.contains_key(&cluster_id) {
            let path = cluster_fastq_path(&clusters_dir, cluster_id);
            let file = File::create(&path)
                .with_context(|| format!("unable to create {}", path.display()))?;
            writers.insert(cluster_id, BufWriter::new(file));
        }
        let writer = writers.get_mut(&cluster_id).unwrap();
        writer.write_all(b"@")?;
        writer.write_all(read_id.as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&seq.sequence)?;
        writer.write_all(b"\n+\n")?;
        let qual = encode_fastq_quality(&seq.qualities, seq.sequence.len(), &mut quality_buffer);
        writer.write_all(qual)?;
        writer.write_all(b"\n")?;
        if (idx + 1) % 1_000 == 0 {
            eprintln!("split_reads: wrote {} reads to cluster FASTQs...", idx + 1);
        }
    }

    Ok(())
}

fn cluster_fastq_path(dir: &Path, cluster_id: i32) -> PathBuf {
    if cluster_id >= 0 {
        dir.join(format!("cluster_{cluster_id}.fastq"))
    } else {
        dir.join("cluster_noise.fastq")
    }
}

fn write_cluster_bams(
    bam_path: &Path,
    clusters_dir: PathBuf,
    read_ids: &[String],
    assignments: &[i32],
) -> Result<()> {
    fs::create_dir_all(&clusters_dir)
        .with_context(|| format!("failed to create {}", clusters_dir.display()))?;

    let mut cluster_lookup: HashMap<&str, i32> = HashMap::new();
    for (read_id, cluster_id) in read_ids.iter().zip(assignments.iter()) {
        cluster_lookup.insert(read_id.as_str(), *cluster_id);
    }

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_path)
        .with_context(|| format!("failed to open {}", bam_path.display()))?;
    let header = reader
        .read_header()
        .with_context(|| format!("failed to read BAM header from {}", bam_path.display()))?;

    let mut record = Record::default();
    let mut writers: HashMap<i32, bam::io::Writer<bgzf::io::Writer<File>>> = HashMap::new();

    loop {
        let bytes = reader
            .read_record(&mut record)
            .context("failed to read BAM record")?;
        if bytes == 0 {
            break;
        }

        let Some(name) = record.name().map(|n| n.to_str_lossy().into_owned()) else {
            continue;
        };
        let Some(cluster_id) = cluster_lookup.get(name.as_str()) else {
            continue;
        };

        let writer = match writers.entry(*cluster_id) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                let path = cluster_bam_path(&clusters_dir, *cluster_id);
                let mut writer = bam::io::writer::Builder::default()
                    .build_from_path(&path)
                    .with_context(|| format!("failed to create {}", path.display()))?;
                writer
                    .write_header(&header)
                    .with_context(|| format!("failed to write BAM header to {}", path.display()))?;
                entry.insert(writer)
            }
        };

        writer
            .write_record(&header, &record)
            .context("failed to write cluster BAM record")?;
    }

    Ok(())
}

fn cluster_bam_path(dir: &Path, cluster_id: i32) -> PathBuf {
    if cluster_id >= 0 {
        dir.join(format!("cluster_{cluster_id}.bam"))
    } else {
        dir.join("cluster_noise.bam")
    }
}

fn auto_min_cluster_size(read_count: usize) -> usize {
    let minimum = 25;
    if read_count == 0 {
        return minimum;
    }
    ((read_count as f64 / 50.0).ceil() as usize).max(minimum)
}

fn auto_hdbscan_min_samples(read_count: usize) -> usize {
    let minimum = 5;
    if read_count == 0 {
        return minimum;
    }
    ((read_count as f64 / 50.0).ceil() as usize).max(minimum)
}

fn encode_fastq_quality<'a>(qualities: &[u8], len: usize, buffer: &'a mut Vec<u8>) -> &'a [u8] {
    buffer.clear();
    if qualities.is_empty() {
        buffer.resize(len, b'!');
        return buffer.as_slice();
    }
    buffer.reserve(len);
    for q in qualities.iter().copied().take(len) {
        buffer.push(q.saturating_add(33));
    }
    while buffer.len() < len {
        buffer.push(b'!');
    }
    buffer.as_slice()
}

fn retain_dense_features(
    rows: Vec<Vec<Option<f64>>>,
    motif_names: Vec<String>,
    max_missing_fraction: f64,
) -> Result<(Vec<Vec<Option<f64>>>, Vec<String>, usize)> {
    if rows.is_empty() || motif_names.is_empty() {
        return Ok((rows, motif_names, 0));
    }

    let row_count = rows.len();
    let mut keep = vec![true; motif_names.len()];
    for col_idx in 0..motif_names.len() {
        let missing = rows.iter().filter(|row| row[col_idx].is_none()).count();
        if row_count > 0 && (missing as f64) / row_count as f64 > max_missing_fraction {
            keep[col_idx] = false;
        }
    }

    let dropped = keep.iter().filter(|flag| !**flag).count();
    if dropped == 0 {
        return Ok((rows, motif_names, 0));
    }

    let mut filtered_rows = Vec::with_capacity(rows.len());
    for row in rows {
        let mut filtered = Vec::with_capacity(motif_names.len() - dropped);
        for (idx, value) in row.into_iter().enumerate() {
            if keep[idx] {
                filtered.push(value);
            }
        }
        filtered_rows.push(filtered);
    }

    let filtered_names = motif_names
        .into_iter()
        .enumerate()
        .filter_map(|(idx, name)| keep[idx].then_some(name))
        .collect();

    Ok((filtered_rows, filtered_names, dropped))
}

fn drop_reads_without_features(
    read_ids: Vec<String>,
    features: Vec<Vec<Option<f64>>>,
    mut sequences: HashMap<String, SequenceInfo>,
) -> (
    Vec<String>,
    Vec<Vec<Option<f64>>>,
    HashMap<String, SequenceInfo>,
) {
    let mut keep_ids = Vec::new();
    let mut keep_features = Vec::new();
    let mut keep_set = HashSet::new();

    for (read_id, row) in read_ids.into_iter().zip(features.into_iter()) {
        if row.iter().any(|v| v.is_some()) {
            keep_set.insert(read_id.clone());
            keep_ids.push(read_id);
            keep_features.push(row);
        }
    }

    sequences.retain(|id, _| keep_set.contains(id));
    (keep_ids, keep_features, sequences)
}

fn impute_features(rows: &[Vec<Option<f64>>]) -> Vec<Vec<f64>> {
    if rows.is_empty() {
        return Vec::new();
    }

    let row_count = rows.len();
    let motif_count = rows[0].len();
    let mut imputed = vec![vec![0.0f64; motif_count]; row_count];
    let mut row_entries: Vec<Vec<(usize, f64)>> = Vec::with_capacity(row_count);
    let mut column_entries: Vec<Vec<(usize, f64)>> = vec![Vec::new(); motif_count];
    let mut column_sums = vec![0.0f64; motif_count];
    let mut column_counts = vec![0usize; motif_count];

    for (row_idx, row) in rows.iter().enumerate() {
        let mut entries = Vec::new();
        for (col_idx, value) in row.iter().enumerate() {
            match value {
                Some(v) => {
                    imputed[row_idx][col_idx] = *v;
                    entries.push((col_idx, *v));
                    column_entries[col_idx].push((row_idx, *v));
                    column_sums[col_idx] += *v;
                    column_counts[col_idx] += 1;
                }
                None => imputed[row_idx][col_idx] = f64::NAN,
            }
        }
        row_entries.push(entries);
    }

    let mut column_means = vec![0.5f64; motif_count];
    for idx in 0..motif_count {
        if column_counts[idx] > 0 {
            column_means[idx] = column_sums[idx] / column_counts[idx] as f64;
        }
    }

    for row_idx in 0..row_count {
        for col_idx in 0..motif_count {
            if imputed[row_idx][col_idx].is_nan() {
                let value = nearest_value(
                    row_idx,
                    col_idx,
                    &row_entries,
                    &column_entries[col_idx],
                    column_means[col_idx],
                );
                imputed[row_idx][col_idx] = value;
            }
        }
    }

    imputed
}

fn nearest_value(
    target_row: usize,
    column: usize,
    rows: &[Vec<(usize, f64)>],
    column_rows: &[(usize, f64)],
    fallback: f64,
) -> f64 {
    if column_rows.is_empty() {
        return fallback;
    }

    let mut best_dist = f64::INFINITY;
    let mut best_value = None;
    let target = &rows[target_row];

    for &(candidate_idx, candidate_value) in column_rows {
        if candidate_idx == target_row {
            continue;
        }
        let candidate = &rows[candidate_idx];
        if let Some(dist) = row_distance_without_column(target, candidate, column) {
            if dist < best_dist {
                best_dist = dist;
                best_value = Some(candidate_value);
            }
        }
    }

    best_value.unwrap_or(fallback)
}

fn row_distance_without_column(
    a: &[(usize, f64)],
    b: &[(usize, f64)],
    skip_col: usize,
) -> Option<f64> {
    let mut ia = 0usize;
    let mut ib = 0usize;
    let mut sum = 0.0;
    let mut count = 0usize;

    while ia < a.len() && ib < b.len() {
        let (col_a, val_a) = a[ia];
        let (col_b, val_b) = b[ib];
        match col_a.cmp(&col_b) {
            std::cmp::Ordering::Equal => {
                if col_a != skip_col {
                    let diff = val_a - val_b;
                    sum += diff * diff;
                    count += 1;
                }
                ia += 1;
                ib += 1;
            }
            std::cmp::Ordering::Less => ia += 1,
            std::cmp::Ordering::Greater => ib += 1,
        }
    }

    if count == 0 {
        None
    } else {
        Some((sum / count as f64).sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::ClusterAlgorithm;

    #[test]
    fn impute_features_fills_missing_values() {
        let features = vec![
            vec![Some(0.1), Some(0.2)],
            vec![Some(0.3), None],
            vec![Some(0.2), Some(0.25)],
        ];

        let imputed = impute_features(&features);
        assert_eq!(imputed.len(), 3);
        assert_eq!(imputed[0], vec![0.1, 0.2]);
        assert_eq!(imputed[2], vec![0.2, 0.25]);
        assert!((imputed[1][1] - 0.25).abs() < 1e-6);
    }

    #[test]
    fn cluster_reads_assigns_clusters() {
        let data = vec![
            vec![0.0f64, 0.0],
            vec![0.1, 0.1],
            vec![5.0, 5.0],
            vec![5.1, 5.1],
        ];
        let hdbscan = cluster_reads(&data, ClusterAlgorithm::Hdbscan, 2, 2)
            .expect("hdbscan clustering failed");
        assert_eq!(hdbscan.len(), 4);
        assert!(hdbscan[0] != hdbscan[2]);

        let gmm = cluster_reads(&data, ClusterAlgorithm::Gmm, 2, 2).expect("gmm clustering failed");
        assert_eq!(gmm.len(), 4);
        assert!(gmm[0] != gmm[2]);

        let agg = cluster_reads(&data, ClusterAlgorithm::Agglomerative, 2, 2)
            .expect("agglomerative clustering failed");
        assert_eq!(agg.len(), 4);
        assert!(agg[0] != agg[2]);
    }

    #[test]
    fn retain_dense_features_filters_sparse_columns() {
        let rows = vec![
            vec![Some(0.1), None, Some(0.2)],
            vec![Some(0.2), None, None],
            vec![Some(0.3), None, Some(0.4)],
            vec![Some(0.4), None, None],
        ];
        let motifs = vec!["A".into(), "B".into(), "C".into()];
        let (filtered, names, dropped) = retain_dense_features(rows, motifs, 0.5).unwrap();
        assert_eq!(names, vec!["A", "C"]);
        assert_eq!(filtered[0].len(), 2);
        assert_eq!(dropped, 1);
    }
}
