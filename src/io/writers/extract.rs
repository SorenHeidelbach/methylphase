use crate::core::MotifQuery;
use crate::features::extract::{
    AggregateRecord, ExtractionSink, FastqRecord, MotifSummaryRecord, PerReadRecord, QUANTILE_PCTS,
};
use anyhow::{Context, Result};
use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

pub struct ExtractionWriter {
    per_read: BufWriter<File>,
    aggregate: BufWriter<File>,
    motif_summary: Option<BufWriter<File>>,
    fastq_dir: Option<PathBuf>,
    fastq_writers: HashMap<String, BufWriter<File>>,
    per_read_count: usize,
    aggregate_count: usize,
}

impl ExtractionWriter {
    pub fn new<P, Q, R, S>(
        motifs: &[MotifQuery],
        per_read_path: P,
        aggregate_path: Q,
        motif_summary_path: Option<R>,
        fastq_dir: Option<S>,
    ) -> Result<Self>
    where
        P: AsRef<Path>,
        Q: AsRef<Path>,
        R: AsRef<Path>,
        S: AsRef<Path>,
    {
        let per_read_file = File::create(per_read_path.as_ref())
            .with_context(|| format!("unable to create {}", per_read_path.as_ref().display()))?;
        let aggregate_file = File::create(aggregate_path.as_ref())
            .with_context(|| format!("unable to create {}", aggregate_path.as_ref().display()))?;

        let motif_summary = if let Some(path) = motif_summary_path {
            let file = File::create(path.as_ref())
                .with_context(|| format!("unable to create {}", path.as_ref().display()))?;
            let mut writer = BufWriter::new(file);
            write!(writer, "read_id\tcontig")?;
            for label in motifs.iter().map(|m| sanitize_label(m.raw())) {
                write!(writer, "\t{label}")?;
            }
            writeln!(writer)?;
            Some(writer)
        } else {
            None
        };

        let fastq_dir = if let Some(dir) = fastq_dir {
            let path = dir.as_ref();
            fs::create_dir_all(path)
                .with_context(|| format!("failed to create {}", path.display()))?;
            Some(path.to_path_buf())
        } else {
            None
        };

        let mut writer = Self {
            per_read: BufWriter::new(per_read_file),
            aggregate: BufWriter::new(aggregate_file),
            motif_summary,
            fastq_dir,
            fastq_writers: HashMap::new(),
            per_read_count: 0,
            aggregate_count: 0,
        };

        writer.write_headers()?;
        Ok(writer)
    }

    fn write_headers(&mut self) -> Result<()> {
        writeln!(
            self.per_read,
            "read_id\tcontig\tstrand\tmotif\tmotif_start\tmotif_position\tread_position\tref_position\tprobability\tmethylated\tmod_label"
        )?;
        write!(
            self.aggregate,
            "read_id\tmotif\tcall_count\tmean_probability"
        )?;
        for pct in QUANTILE_PCTS {
            write!(self.aggregate, "\tquantile_{pct}")?;
        }
        writeln!(self.aggregate)?;
        Ok(())
    }

    fn fastq_writer(&mut self, key: &str) -> Result<&mut BufWriter<File>> {
        if !self.fastq_writers.contains_key(key) {
            let dir = self.fastq_dir.as_ref().expect("fastq dir exists");
            let path = dir.join(format!("{key}.fastq"));
            let file = File::create(&path)
                .with_context(|| format!("unable to create {}", path.display()))?;
            self.fastq_writers
                .insert(key.to_string(), BufWriter::new(file));
        }

        Ok(self.fastq_writers.get_mut(key).unwrap())
    }
}

impl ExtractionSink for ExtractionWriter {
    fn record_per_read(&mut self, record: PerReadRecord<'_>) -> Result<()> {
        let contig = record.contig.unwrap_or(".");
        let strand = record
            .strand
            .map(|s| s.to_string())
            .unwrap_or_else(|| ".".to_string());
        let ref_position = record
            .reference_position
            .map(|p| p.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let probability = record
            .probability
            .map(|p| format!("{:.4}", p))
            .unwrap_or_else(|| "0.0000".to_string());
        let methylated = if record.methylated { "1" } else { "0" };

        writeln!(
            self.per_read,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.read_id,
            contig,
            strand,
            record.motif_name,
            record.motif_start,
            record.motif_position,
            record.read_position,
            ref_position,
            probability,
            methylated,
            record.mod_label,
        )?;
        self.per_read_count += 1;
        if self.per_read_count % 10_000 == 0 {
            self.per_read.flush()?;
        }

        Ok(())
    }

    fn record_aggregate(&mut self, record: AggregateRecord<'_>) -> Result<()> {
        let probability = record
            .mean_probability
            .map(|p| format!("{:.4}", p))
            .unwrap_or_else(|| "NA".to_string());
        let mut quantiles = Vec::new();
        if let Some(values) = record.quantiles.as_ref() {
            quantiles.extend(values.iter().map(|v| format!("{:.4}", v)));
        } else {
            quantiles.resize(QUANTILE_PCTS.len(), "NA".to_string());
        }

        write!(
            self.aggregate,
            "{}\t{}\t{}\t{}",
            record.read_id, record.motif_name, record.call_count, probability
        )?;
        for value in quantiles {
            write!(self.aggregate, "\t{}", value)?;
        }
        writeln!(self.aggregate)?;
        self.aggregate_count += 1;
        if self.aggregate_count % 10_000 == 0 {
            self.aggregate.flush()?;
        }

        Ok(())
    }

    fn record_motif_summary(&mut self, record: MotifSummaryRecord<'_>) -> Result<()> {
        if let Some(writer) = self.motif_summary.as_mut() {
            write!(
                writer,
                "{}\t{}",
                record.read_id,
                record.contig.unwrap_or(".")
            )?;
            for hit in record.hits {
                write!(writer, "\t{}", if *hit { 1 } else { 0 })?;
            }
            writeln!(writer)?;
        }
        Ok(())
    }

    fn record_fastq(&mut self, record: FastqRecord<'_>) -> Result<()> {
        let Some(_) = self.fastq_dir else {
            return Ok(());
        };

        let writer = self.fastq_writer(record.file_key)?;
        let seq = std::str::from_utf8(record.sequence).unwrap_or("");
        let qual = encode_quality(record.qualities, record.sequence.len());

        writeln!(writer, "@{}|{}", record.read_id, record.combo_label)?;
        writeln!(writer, "{seq}")?;
        writeln!(writer, "+")?;
        writeln!(writer, "{qual}")?;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.per_read.flush()?;
        self.aggregate.flush()?;
        if let Some(writer) = self.motif_summary.as_mut() {
            writer.flush()?;
        }
        for writer in self.fastq_writers.values_mut() {
            writer.flush()?;
        }
        Ok(())
    }
}

fn encode_quality(qualities: &[u8], len: usize) -> String {
    if qualities.is_empty() {
        return std::iter::repeat('!').take(len).collect::<String>();
    }

    qualities
        .iter()
        .map(|q| (q.saturating_add(33)) as char)
        .collect()
}

fn sanitize_label(raw: &str) -> String {
    raw.chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect()
}
