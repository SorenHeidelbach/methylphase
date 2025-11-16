use crate::core::{MotifQuery, ReadRecord, Strand};
use anyhow::{bail, Result};
use std::collections::HashMap;

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
}

impl<S: ExtractionSink> Extractor<S> {
    pub fn new(motifs: Vec<MotifQuery>, sink: S) -> Result<Self> {
        if motifs.is_empty() {
            bail!("at least one motif must be supplied");
        }

        Ok(Self { motifs, sink })
    }

    pub fn process_read(&mut self, read: &ReadRecord) -> Result<()> {
        let mut per_read_stats: HashMap<&str, MotifStats> = HashMap::new();
        let mut motif_hits = vec![false; self.motifs.len()];

        for (idx, motif) in self.motifs.iter().enumerate() {
            let motif_name = motif.raw();
            let motif_position = motif.mod_position();
            let mod_label = motif.mod_label();
            let motif_offset = motif.mod_position() - 1;
            let matches = motif.matches(read.sequence());

            for start in matches {
                let mod_pos = start + motif_offset;
                let Some(call) = read.modification_at(mod_pos, mod_label) else {
                    continue;
                };

                motif_hits[idx] = true;
                self.sink.record_per_read(PerReadRecord {
                    read_id: &read.id,
                    contig: read.contig_name(),
                    strand: read.strand,
                    motif_name,
                    motif_start: start,
                    motif_position,
                    read_position: mod_pos,
                    reference_position: read.reference_position(mod_pos),
                    probability: call.probability,
                    mod_label: &call.kind.label,
                })?;

                per_read_stats
                    .entry(motif_name)
                    .or_default()
                    .observe(call.probability);
            }
        }

        for (motif_name, stats) in per_read_stats {
            self.sink.record_aggregate(AggregateRecord {
                read_id: &read.id,
                motif_name,
                call_count: stats.count,
                mean_probability: stats.mean_probability(),
                quantiles: stats.quantiles(),
            })?;
        }

        self.sink.record_motif_summary(MotifSummaryRecord {
            read_id: &read.id,
            contig: read.contig_name(),
            hits: &motif_hits,
        })?;

        let (combo_label, file_key) = motif_combo_key(&motif_hits, &self.motifs);
        self.sink.record_fastq(FastqRecord {
            read_id: &read.id,
            combo_label: &combo_label,
            file_key: &file_key,
            sequence: read.sequence(),
            qualities: read.quality_scores(),
        })?;

        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.sink.finish()
    }
}

#[derive(Default)]
struct MotifStats {
    count: usize,
    probability_sum: f64,
    probabilities: Vec<f32>,
}

impl MotifStats {
    fn observe(&mut self, probability: Option<f32>) {
        self.count += 1;
        if let Some(p) = probability {
            self.probability_sum += f64::from(p);
            self.probabilities.push(p);
        }
    }

    fn mean_probability(&self) -> Option<f64> {
        (!self.probabilities.is_empty()).then_some(self.probability_sum / self.probabilities.len() as f64)
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
    pub mod_label: &'a str,
}

pub struct AggregateRecord<'a> {
    pub read_id: &'a str,
    pub motif_name: &'a str,
    pub call_count: usize,
    pub mean_probability: Option<f64>,
    pub quantiles: Option<[f64; QUANTILE_COUNT]>,
}

pub struct MotifSummaryRecord<'a> {
    pub read_id: &'a str,
    pub contig: Option<&'a str>,
    pub hits: &'a [bool],
}

pub struct FastqRecord<'a> {
    pub read_id: &'a str,
    pub combo_label: &'a str,
    pub file_key: &'a str,
    pub sequence: &'a [u8],
    pub qualities: &'a [u8],
}

fn motif_combo_key(hits: &[bool], motifs: &[MotifQuery]) -> (String, String) {
    let names: Vec<&str> = hits
        .iter()
        .enumerate()
        .filter_map(|(idx, hit)| (*hit).then_some(motifs[idx].raw()))
        .collect();

    let label = if names.is_empty() {
        "NONE".to_string()
    } else {
        names.join("+")
    };

    let file_key = format!("combo_{}", sanitize_label(&label));
    (label, file_key)
}

fn sanitize_label(raw: &str) -> String {
    raw.chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect()
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
