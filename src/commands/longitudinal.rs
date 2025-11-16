use crate::{
    commands::shared,
    core::{MotifQuery, ReadRecord},
    features::longitudinal::{BlockResult, LongitudinalAnalyzer},
    io::{BamReader, BamRecordDecoder},
};
use anyhow::{bail, Context, Result};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

pub fn run(
    bam: PathBuf,
    motifs: Vec<String>,
    motif_file: Option<PathBuf>,
    output: PathBuf,
    block_size: usize,
    methylation_threshold: f32,
    contigs: Vec<String>,
) -> Result<()> {
    if block_size == 0 {
        bail!("block-size must be greater than zero");
    }

    let motif_queries = shared::load_motif_queries(motifs, motif_file)?;
    let block_size_bp = i64::try_from(block_size).context("block size exceeds i64 range")?;
    let mut analyzer = LongitudinalAnalyzer::new(block_size_bp);

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
        process_read(&read, &motif_queries, methylation_threshold, &mut analyzer);
        Ok(())
    })?;

    decoder.report();
    let results = analyzer.finish();
    write_results(&results, &output)?;

    Ok(())
}

fn process_read(
    read: &ReadRecord,
    motifs: &[MotifQuery],
    threshold: f32,
    analyzer: &mut LongitudinalAnalyzer,
) {
    let Some(contig) = read.contig_name() else {
        return;
    };

    for motif in motifs {
        let motif_offset = motif.mod_position() - 1;
        let matches = motif.matches(read.sequence());
        for start in matches {
            let mod_pos = start + motif_offset;
            let Some(ref_pos) = read.reference_position(mod_pos) else {
                continue;
            };
            let state = read
                .modification_at(mod_pos, motif.mod_label())
                .map(|call| call.probability.map(|p| p >= threshold).unwrap_or(true))
                .unwrap_or(false);
            analyzer.record_call(contig, ref_pos, motif.raw(), &read.id, state);
        }
    }
}

fn write_results(results: &[BlockResult], output: &PathBuf) -> Result<()> {
    let file =
        File::create(output).with_context(|| format!("unable to create {}", output.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "contig\tblock_start\tblock_end\tsite_count\tpattern\tread_count\tsite_positions\t\
         total_reads\tinformative_reads\tuninformative_reads\tfully_methylated_reads\tfully_unmethylated_reads\tmixed_reads\t\
         mean_methylation\tmethylation_variance\tpolarization"
    )?;

    for block in results {
        let site_positions = block
            .sites
            .iter()
            .map(|site| format!("{}:{}", site.position, site.motif))
            .collect::<Vec<_>>()
            .join(",");
        for pattern in &block.patterns {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                block.contig,
                block.block_start,
                block.block_end,
                block.sites.len(),
                pattern.pattern,
                pattern.read_count,
                site_positions,
                block.statistics.total_reads,
                block.statistics.informative_reads,
                block.statistics.uninformative_reads,
                block.statistics.fully_methylated_reads,
                block.statistics.fully_unmethylated_reads,
                block.statistics.mixed_reads,
                block
                    .statistics
                    .mean_methylation
                    .map(|v| format!("{:.4}", v))
                    .unwrap_or_else(|| "NA".to_string()),
                block
                    .statistics
                    .methylation_variance
                    .map(|v| format!("{:.4}", v))
                    .unwrap_or_else(|| "NA".to_string()),
                block
                    .statistics
                    .polarization
                    .map(|v| format!("{:.4}", v))
                    .unwrap_or_else(|| "NA".to_string()),
            )?;
        }
    }

    Ok(())
}
