use crate::{
    commands::shared,
    features::extract::Extractor,
    io::{writers::ExtractionWriter, BamReader, BamRecordDecoder, SequenceCache},
};
use anyhow::{Context, Result};
use std::{
    fs,
    path::{Path, PathBuf},
    sync::Arc,
};

#[allow(clippy::too_many_arguments)]
pub fn run(
    bam: PathBuf,
    motifs: Vec<String>,
    motif_file: Option<PathBuf>,
    motif_summary_tsv: Option<PathBuf>,
    fastq_dir: Option<PathBuf>,
    sequence_fallback: Option<PathBuf>,
    sequence_index: Option<PathBuf>,
    output_dir: PathBuf,
    contigs: Vec<String>,
) -> Result<()> {
    let motif_queries = shared::load_motif_queries(motifs, motif_file)?;
    let motif_count = motif_queries.len();

    fs::create_dir_all(&output_dir)
        .with_context(|| format!("failed to create {}", output_dir.display()))?;
    let outputs = resolve_outputs(&output_dir, motif_summary_tsv, fastq_dir)?;
    eprintln!(
        "extract: processing {} motifs; writing outputs to {}",
        motif_count,
        output_dir.display()
    );
    let writer = ExtractionWriter::new(
        &motif_queries,
        outputs.per_read,
        outputs.aggregate,
        outputs.motif_summary,
        outputs.fastq_dir,
    )?;
    let mut extractor = Extractor::new(motif_queries, writer)?;

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

    let mut processed_reads = 0usize;
    reader.visit_records(contig_filter.as_ref().map(|c| c.as_slice()), |record| {
        let read = decoder.decode(record)?;
        extractor.process_read(&read)?;
        processed_reads += 1;
        if processed_reads % 1_000 == 0 {
            eprintln!("extract: processed {} reads...", processed_reads);
        }
        Ok(())
    })?;

    decoder.report();
    extractor.finish()?;
    eprintln!("extract: finished, processed {} reads", processed_reads);
    Ok(())
}

struct ExtractOutputs {
    per_read: PathBuf,
    aggregate: PathBuf,
    motif_summary: Option<PathBuf>,
    fastq_dir: Option<PathBuf>,
}

fn resolve_outputs(
    output_dir: &Path,
    motif_summary: Option<PathBuf>,
    fastq_dir: Option<PathBuf>,
) -> Result<ExtractOutputs> {
    let per_read_path = output_dir.join("per_read.tsv");
    let aggregate_path = output_dir.join("aggregate.tsv");
    let motif_summary_path = motif_summary.map(|path| join_relative(Some(output_dir), path));
    let fastq_dir_path = fastq_dir.map(|path| join_relative(Some(output_dir), path));
    Ok(ExtractOutputs {
        per_read: per_read_path,
        aggregate: aggregate_path,
        motif_summary: motif_summary_path,
        fastq_dir: fastq_dir_path,
    })
}

fn join_relative(base: Option<&Path>, path: PathBuf) -> PathBuf {
    if let Some(dir) = base {
        if path.is_relative() {
            return dir.join(path);
        }
    }
    path
}
