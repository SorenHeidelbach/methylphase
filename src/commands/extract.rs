use crate::{
    commands::shared,
    features::extract::Extractor,
    io::writers::ExtractionWriter,
    io::{BamReader, BamRecordDecoder},
};
use anyhow::{Context, Result};
use std::path::PathBuf;

#[allow(clippy::too_many_arguments)]
pub fn run(
    bam: PathBuf,
    motifs: Vec<String>,
    motif_file: Option<PathBuf>,
    motif_summary_tsv: Option<PathBuf>,
    fastq_dir: Option<PathBuf>,
    per_read_tsv: PathBuf,
    aggregate_tsv: PathBuf,
    contigs: Vec<String>,
) -> Result<()> {
    let motif_queries = shared::load_motif_queries(motifs, motif_file)?;

    let writer = ExtractionWriter::new(
        &motif_queries,
        per_read_tsv,
        aggregate_tsv,
        motif_summary_tsv,
        fastq_dir,
    )?;
    let mut extractor = Extractor::new(motif_queries, writer)?;

    let mut reader =
        BamReader::open(&bam).with_context(|| format!("failed to open BAM {}", bam.display()))?;
    let decoder = BamRecordDecoder::new(reader.header());

    let contig_filter = if contigs.is_empty() {
        None
    } else {
        Some(contigs)
    };

    reader.visit_records(contig_filter.as_ref().map(|c| c.as_slice()), |record| {
        let read = decoder.decode(record)?;
        extractor.process_read(&read)
    })?;

    decoder.report();
    extractor.finish()?;
    Ok(())
}
