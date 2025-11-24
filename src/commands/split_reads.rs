use crate::{
    commands::shared,
    core::{MotifQuery, ReadRecord},
    features::extract::{
        self, AggregateRecord, ExtractionSink, FastqRecord, MotifSummaryRecord, PerReadRecord,
    },
    io::{writers::ExtractionWriter, BamReader, BamRecordDecoder, SequenceCache},
};
use anyhow::{anyhow, bail, Context, Result};
use bstr::ByteSlice;
use crossbeam_channel::{bounded, Receiver, Sender};
use hdbscan::{Hdbscan, HdbscanHyperParams};
use noodles_bam::{self as bam, Record};
use noodles_bgzf as bgzf;
use std::{
    collections::{hash_map::Entry, HashMap},
    fs::{self, File},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
    thread,
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
    contigs: Vec<String>,
) -> Result<()> {
    if min_cluster_size == 0 {
        bail!("min-cluster-size must be greater than zero");
    }
    let motif_queries = shared::load_motif_queries(motifs, motif_file)?;
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
        let count = process_reads_parallel(
            &mut reader,
            &decoder,
            contig_filter.as_ref().map(|c| c.as_slice()),
            motif_queries.clone(),
            shared.clone(),
            threads,
        )?;
        let sink = shared.into_inner()?;
        (sink, count)
    };
    decoder.report();
    eprintln!(
        "split_reads: finished extraction phase ({} reads); starting clustering",
        processed_reads
    );
    let (reads, raw_features, sequences, motif_names) = sink.into_data();
    if reads.is_empty() {
        bail!("no reads matched motifs; nothing to cluster");
    }
    eprintln!(
        "Imputing feature matrix for {} reads and {} motifs...",
        reads.len(), motif_names.len()
    );
    let imputed_features = impute_features(&raw_features);

    let effective_min_cluster_size = auto_min_cluster_size(reads.len());
    let effective_min_samples = auto_hdbscan_min_samples(reads.len());
    eprintln!(
        "split_reads: clustering {} reads (min_cluster_size={}, min_samples={})",
        reads.len(),
        effective_min_cluster_size,
        effective_min_samples
    );
    let assignments =
        cluster_reads(&imputed_features, effective_min_cluster_size, effective_min_samples)?;
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

struct SharedSink<S> {
    inner: Arc<Mutex<S>>,
}

impl<S> Clone for SharedSink<S> {
    fn clone(&self) -> Self {
        Self {
            inner: Arc::clone(&self.inner),
        }
    }
}

impl<S> SharedSink<S> {
    fn new(inner: S) -> Self {
        Self {
            inner: Arc::new(Mutex::new(inner)),
        }
    }

    fn into_inner(self) -> Result<S> {
        match Arc::try_unwrap(self.inner) {
            Ok(mutex) => mutex
                .into_inner()
                .map_err(|e| anyhow!("failed to unlock sink: {e}")),
            Err(_) => Err(anyhow!("sink still has outstanding references at shutdown")),
        }
    }
}

impl<S: ExtractionSink> ExtractionSink for SharedSink<S> {
    fn record_per_read(&mut self, record: PerReadRecord<'_>) -> Result<()> {
        let mut guard = self.inner.lock().map_err(|e| anyhow!(e.to_string()))?;
        guard.record_per_read(record)
    }

    fn record_aggregate(&mut self, record: AggregateRecord<'_>) -> Result<()> {
        let mut guard = self.inner.lock().map_err(|e| anyhow!(e.to_string()))?;
        guard.record_aggregate(record)
    }

    fn record_motif_summary(&mut self, record: MotifSummaryRecord<'_>) -> Result<()> {
        let mut guard = self.inner.lock().map_err(|e| anyhow!(e.to_string()))?;
        guard.record_motif_summary(record)
    }

    fn record_fastq(&mut self, record: FastqRecord<'_>) -> Result<()> {
        let mut guard = self.inner.lock().map_err(|e| anyhow!(e.to_string()))?;
        guard.record_fastq(record)
    }

    fn finish(&mut self) -> Result<()> {
        let mut guard = self.inner.lock().map_err(|e| anyhow!(e.to_string()))?;
        guard.finish()
    }
}

struct ClusterSink<W: ExtractionSink> {
    inner: W,
    motif_names: Vec<String>,
    motif_index: HashMap<String, usize>,
    read_features: HashMap<String, Vec<Option<f64>>>,
    sequences: HashMap<String, SequenceInfo>,
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
        let mut read_ids: Vec<String> = self.read_features.keys().cloned().collect();
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
        let qual = encode_fastq_quality(&seq.qualities, seq.sequence.len());
        writer.write_all(qual.as_bytes())?;
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
    ((read_count as f64).sqrt().ceil() as usize).max(minimum)
}

fn auto_hdbscan_min_samples(read_count: usize) -> usize {
    let minimum = 5;
    if read_count == 0 {
        return minimum;
    }
    ((read_count as f64 / 100.0).ceil() as usize).max(minimum)
}

fn process_reads_parallel<S>(
    reader: &mut BamReader,
    decoder: &BamRecordDecoder,
    contigs: Option<&[String]>,
    motifs: Vec<MotifQuery>,
    shared_sink: SharedSink<S>,
    threads: usize,
) -> Result<usize>
where
    S: ExtractionSink + Send + 'static,
{
    let (tx, rx): (Sender<ReadRecord>, Receiver<ReadRecord>) = bounded(threads * 4);
    let processed = Arc::new(AtomicUsize::new(0));
    let mut handles = Vec::new();
    for _ in 0..threads {
        let rx = rx.clone();
        let motifs = motifs.clone();
        let mut sink = shared_sink.clone();
        let processed = processed.clone();
        handles.push(thread::spawn(move || -> Result<usize> {
            let mut local = 0usize;
            while let Ok(read) = rx.recv() {
                extract::process_read_with_motifs(&read, &motifs, &mut sink)?;
                local += 1;
                let total = processed.fetch_add(1, Ordering::Relaxed) + 1;
                if total % 1_000 == 0 {
                    eprintln!("split_reads: processed {} reads...", total);
                }
            }
            Ok(local)
        }));
    }
    drop(rx);

    reader.visit_records(contigs, |record| {
        let read = decoder.decode(record)?;
        tx.send(read)
            .map_err(|e| anyhow!("failed to dispatch read: {e}"))?;
        Ok(())
    })?;
    drop(tx);

    let mut total = 0usize;
    for handle in handles {
        total += handle
            .join()
            .map_err(|_| anyhow!("worker thread panicked"))??;
    }
    Ok(total)
}

fn encode_fastq_quality(qualities: &[u8], len: usize) -> String {
    if qualities.is_empty() {
        return std::iter::repeat('!').take(len).collect();
    }
    qualities
        .iter()
        .map(|q| (q.saturating_add(33)) as char)
        .collect()
}

fn impute_features(rows: &[Vec<Option<f64>>]) -> Vec<Vec<f64>> {
    if rows.is_empty() {
        return Vec::new();
    }

    let motif_count = rows[0].len();
    let mut column_means = vec![0.0; motif_count];
    let mut column_counts = vec![0usize; motif_count];
    for row in rows {
        for (idx, value) in row.iter().enumerate() {
            if let Some(v) = value {
                column_means[idx] += *v;
                column_counts[idx] += 1;
            }
        }
    }
    for idx in 0..motif_count {
        if column_counts[idx] > 0 {
            column_means[idx] /= column_counts[idx] as f64;
        } else {
            column_means[idx] = 0.5;
        }
    }

    rows.iter()
        .enumerate()
        .map(|(row_idx, row)| {
            row.iter()
                .enumerate()
                .map(|(col_idx, value)| match value {
                    Some(v) => *v,
                    None => nearest_value(rows, row_idx, col_idx).unwrap_or(column_means[col_idx]),
                })
                .collect::<Vec<f64>>()
        })
        .collect()
}

fn nearest_value(rows: &[Vec<Option<f64>>], target_row: usize, column: usize) -> Option<f64> {
    let mut best_dist = f64::INFINITY;
    let mut best_value = None;
    for (idx, row) in rows.iter().enumerate() {
        if idx == target_row {
            continue;
        }
        if let Some(value) = row[column] {
            if let Some(dist) = row_distance_without_column(&rows[target_row], row, column) {
                if dist < best_dist {
                    best_dist = dist;
                    best_value = Some(value);
                }
            }
        }
    }
    best_value
}

fn row_distance_without_column(
    a: &[Option<f64>],
    b: &[Option<f64>],
    skip_col: usize,
) -> Option<f64> {
    let mut sum = 0.0;
    let mut count = 0usize;
    for (idx, (av, bv)) in a.iter().zip(b.iter()).enumerate() {
        if idx == skip_col {
            continue;
        }
        // eucledian distance
        if let (Some(x), Some(y)) = (av, bv) {
            let diff = x - y;
            sum += diff * diff;
            count += 1;
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
        let labels = cluster_reads(&data, 2, 2).expect("clustering failed");
        assert_eq!(labels.len(), 4);
        assert!(labels[0] != labels[2]);
    }


}
