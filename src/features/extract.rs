use crate::core::{MotifQuery, ReadRecord, Strand};
use anyhow::{bail, Result};

pub const QUANTILE_COUNT: usize = 11;
pub const QUANTILE_PCTS: [u8; QUANTILE_COUNT] = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];

/// Trait implemented by sinks that persist extraction results.
pub trait ExtractionSink {
    fn record_per_read(&mut self, _record: PerReadRecord<'_>) -> Result<()> {
        Ok(())
    }

    fn record_aggregate(&mut self, _record: AggregateRecord<'_>) -> Result<()> {
        Ok(())
    }

    fn record_motif_summary(&mut self, _record: MotifSummaryRecord<'_>) -> Result<()> {
        Ok(())
    }

    fn record_fastq(&mut self, _record: FastqRecord<'_>) -> Result<()> {
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        Ok(())
    }
}

pub struct Extractor<S: ExtractionSink> {
    motifs: Vec<MotifQuery>,
    sink: S,
    scratch: MotifScratch,
}

impl<S: ExtractionSink> Extractor<S> {
    pub fn new(motifs: Vec<MotifQuery>, sink: S) -> Result<Self> {
        if motifs.is_empty() {
            bail!("at least one motif must be supplied");
        }

        let scratch = MotifScratch::new(motifs.len());
        Ok(Self {
            motifs,
            sink,
            scratch,
        })
    }

    pub fn process_read(&mut self, read: &ReadRecord) -> Result<()> {
        process_read_with_motifs(read, &self.motifs, &mut self.sink, &mut self.scratch)
    }

    pub fn finish(mut self) -> Result<()> {
        self.sink.finish()
    }

    pub fn finish_with_sink(mut self) -> Result<S> {
        self.sink.finish()?;
        Ok(self.sink)
    }
}

pub fn process_read_with_motifs<S: ExtractionSink>(
    read: &ReadRecord,
    motifs: &[MotifQuery],
    sink: &mut S,
    scratch: &mut MotifScratch,
) -> Result<()> {
    scratch.prepare(motifs.len());

    for (idx, motif) in motifs.iter().enumerate() {
        let motif_name = motif.raw();
        let motif_position = motif.mod_position();
        let mod_label = motif.mod_label();
        let motif_offset = motif.mod_position();
        let matches = motif.matches(read.sequence());

        for start in matches {
            let mod_pos = start + motif_offset;
            let Some(call) = read.modification_at(mod_pos, mod_label) else {
                continue;
            };

            if call.methylated {
                scratch.hits[idx] = true;
            }
            sink.record_per_read(PerReadRecord {
                read_id: &read.id,
                contig: read.contig_name(),
                strand: read.strand,
                motif_name,
                motif_start: start,
                motif_position,
                read_position: mod_pos,
                reference_position: read.reference_position(mod_pos),
                probability: call.probability,
                methylated: call.methylated,
                mod_label: &call.kind.label,
            })?;

            scratch.stats[idx].observe(call.probability);
        }
    }

    for (idx, stats) in scratch.stats[..motifs.len()].iter().enumerate() {
        if stats.count == 0 {
            continue;
        }
        sink.record_aggregate(AggregateRecord {
            read_id: &read.id,
            motif_name: motifs[idx].raw(),
            call_count: stats.count,
            mean_probability: stats.mean_probability(),
            quantiles: stats.quantiles(),
        })?;
    }

    sink.record_motif_summary(MotifSummaryRecord {
        read_id: &read.id,
        contig: read.contig_name(),
        hits: &scratch.hits[..motifs.len()],
    })?;

    let (combo_label, file_key) = motif_combo_key(
        &scratch.hits[..motifs.len()],
        motifs,
        &mut scratch.combo_label,
        &mut scratch.file_key,
    );
    sink.record_fastq(FastqRecord {
        read_id: &read.id,
        combo_label,
        file_key,
        sequence: read.sequence(),
        qualities: read.quality_scores(),
    })?;

    Ok(())
}

pub struct MotifScratch {
    hits: Vec<bool>,
    stats: Vec<MotifStats>,
    combo_label: String,
    file_key: String,
}

impl MotifScratch {
    pub fn new(motif_count: usize) -> Self {
        Self {
            hits: vec![false; motif_count],
            stats: vec![MotifStats::default(); motif_count],
            combo_label: String::new(),
            file_key: String::new(),
        }
    }

    fn prepare(&mut self, motif_count: usize) {
        if self.hits.len() < motif_count {
            self.hits.resize(motif_count, false);
        }
        if self.stats.len() < motif_count {
            self.stats.resize_with(motif_count, MotifStats::default);
        }
        for hit in &mut self.hits[..motif_count] {
            *hit = false;
        }
        for stats in &mut self.stats[..motif_count] {
            stats.reset();
        }
    }
}

#[derive(Default, Clone)]
struct MotifStats {
    count: usize,
    probability_sum: f64,
    probabilities: Vec<f32>,
}

impl MotifStats {
    fn reset(&mut self) {
        self.count = 0;
        self.probability_sum = 0.0;
        self.probabilities.clear();
    }

    fn observe(&mut self, probability: Option<f32>) {
        self.count += 1;
        if let Some(p) = probability {
            self.probability_sum += f64::from(p);
            self.probabilities.push(p);
        }
    }

    fn mean_probability(&self) -> Option<f64> {
        (!self.probabilities.is_empty())
            .then_some(self.probability_sum / self.probabilities.len() as f64)
    }

    fn quantiles(&self) -> Option<[f64; QUANTILE_COUNT]> {
        if self.probabilities.is_empty() {
            return None;
        }

        let mut values = self.probabilities.clone();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut result = [0.0; QUANTILE_COUNT];
        for (idx, pct) in QUANTILE_PCTS.iter().enumerate() {
            let fraction = f64::from(*pct) / 100.0;
            result[idx] = percentile(&values, fraction);
        }
        Some(result)
    }
}

pub struct PerReadRecord<'a> {
    pub read_id: &'a str,
    pub contig: Option<&'a str>,
    pub strand: Option<Strand>,
    pub motif_name: &'a str,
    pub motif_start: usize,
    pub motif_position: usize,
    pub read_position: usize,
    pub reference_position: Option<i64>,
    pub probability: Option<f32>,
    pub methylated: bool,
    pub mod_label: &'a str,
}

impl<'a> PerReadRecord<'a> {
    pub fn into_owned(self) -> PerReadRecord<'static> {
        PerReadRecord {
            read_id: Box::leak(self.read_id.to_owned().into_boxed_str()),
            contig: self
                .contig
                .map(|c| Box::leak(c.to_owned().into_boxed_str()) as *const str)
                .map(|ptr| unsafe { &*ptr }),
            strand: self.strand,
            motif_name: Box::leak(self.motif_name.to_owned().into_boxed_str()),
            motif_start: self.motif_start,
            motif_position: self.motif_position,
            read_position: self.read_position,
            reference_position: self.reference_position,
            probability: self.probability,
            methylated: self.methylated,
            mod_label: Box::leak(self.mod_label.to_owned().into_boxed_str()),
        }
    }
}

pub struct AggregateRecord<'a> {
    pub read_id: &'a str,
    pub motif_name: &'a str,
    pub call_count: usize,
    pub mean_probability: Option<f64>,
    pub quantiles: Option<[f64; QUANTILE_COUNT]>,
}

impl<'a> AggregateRecord<'a> {
    pub fn into_owned(self) -> AggregateRecord<'static> {
        AggregateRecord {
            read_id: Box::leak(self.read_id.to_owned().into_boxed_str()),
            motif_name: Box::leak(self.motif_name.to_owned().into_boxed_str()),
            call_count: self.call_count,
            mean_probability: self.mean_probability,
            quantiles: self.quantiles,
        }
    }
}

pub struct MotifSummaryRecord<'a> {
    pub read_id: &'a str,
    pub contig: Option<&'a str>,
    pub hits: &'a [bool],
}

impl<'a> MotifSummaryRecord<'a> {
    pub fn into_owned(self) -> MotifSummaryRecord<'static> {
        MotifSummaryRecord {
            read_id: Box::leak(self.read_id.to_owned().into_boxed_str()),
            contig: self
                .contig
                .map(|c| Box::leak(c.to_owned().into_boxed_str()) as *const str)
                .map(|ptr| unsafe { &*ptr }),
            hits: Box::leak(self.hits.to_vec().into_boxed_slice()),
        }
    }
}

pub struct FastqRecord<'a> {
    pub read_id: &'a str,
    pub combo_label: &'a str,
    pub file_key: &'a str,
    pub sequence: &'a [u8],
    pub qualities: &'a [u8],
}

impl<'a> FastqRecord<'a> {
    pub fn into_owned(self) -> FastqRecord<'static> {
        FastqRecord {
            read_id: Box::leak(self.read_id.to_owned().into_boxed_str()),
            combo_label: Box::leak(self.combo_label.to_owned().into_boxed_str()),
            file_key: Box::leak(self.file_key.to_owned().into_boxed_str()),
            sequence: Box::leak(self.sequence.to_vec().into_boxed_slice()),
            qualities: Box::leak(self.qualities.to_vec().into_boxed_slice()),
        }
    }
}

fn motif_combo_key<'a>(
    hits: &[bool],
    motifs: &[MotifQuery],
    label: &'a mut String,
    file_key: &'a mut String,
) -> (&'a str, &'a str) {
    label.clear();
    file_key.clear();

    let mut any = false;
    for (idx, hit) in hits.iter().enumerate() {
        if *hit {
            if any {
                label.push('+');
            }
            label.push_str(motifs[idx].raw());
            any = true;
        }
    }

    if !any {
        label.push_str("NONE");
    }

    file_key.push_str("combo_");
    for ch in label.chars() {
        if ch.is_ascii_alphanumeric() {
            file_key.push(ch);
        } else {
            file_key.push('_');
        }
    }

    (&*label, &*file_key)
}

fn percentile(values: &[f32], fraction: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }

    let clamped = fraction.clamp(0.0, 1.0);
    if values.len() == 1 {
        return f64::from(values[0]);
    }

    let scaled = clamped * (values.len() - 1) as f64;
    let lower_idx = scaled.floor() as usize;
    let upper_idx = scaled.ceil() as usize;
    if lower_idx == upper_idx {
        f64::from(values[lower_idx])
    } else {
        let lower = f64::from(values[lower_idx]);
        let upper = f64::from(values[upper_idx]);
        let weight = scaled - lower_idx as f64;
        lower + (upper - lower) * weight
    }
}
