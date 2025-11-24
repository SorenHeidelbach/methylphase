pub mod contigs;
pub mod extract;
pub mod impute_bam;
pub mod longitudinal;
pub mod reads;
pub mod shared;
pub mod split_reads;
pub mod vcf;

use crate::cli::{Cli, Command};
use anyhow::Result;

pub fn run(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Contigs { bam } => contigs::run(bam),
        Command::Reads {
            bam,
            contigs,
            limit,
        } => reads::run(bam, contigs, limit),
        Command::Extract {
            bam,
            motifs,
            motif_file,
            motif_summary_tsv,
            fastq_dir,
            sequence_fallback,
            sequence_index,
            output_dir,
            contigs,
        } => extract::run(
            bam,
            motifs,
            motif_file,
            motif_summary_tsv,
            fastq_dir,
            sequence_fallback,
            sequence_index,
            output_dir,
            contigs,
        ),
        Command::Longitudinal {
            bam,
            motifs,
            motif_file,
            block_size,
            methylation_threshold,
            sequence_fallback,
            sequence_index,
            output_dir,
            contigs,
        } => longitudinal::run(
            bam,
            motifs,
            motif_file,
            block_size,
            methylation_threshold,
            sequence_fallback,
            sequence_index,
            output_dir,
            contigs,
        ),
        Command::SplitReads {
            bam,
            motifs,
            motif_file,
            sequence_fallback,
            sequence_index,
            output_dir,
            min_cluster_size,
            min_samples,
            emit_fastq,
            threads,
            contigs,
        } => split_reads::run(
            bam,
            motifs,
            motif_file,
            sequence_fallback,
            sequence_index,
            output_dir,
            min_cluster_size,
            min_samples,
            emit_fastq,
            threads,
            contigs,
        ),
        Command::Vcf {
            bam,
            motifs,
            motif_file,
            sequence_fallback,
            sequence_index,
            output,
            sample_name,
            methylation_threshold,
            plain_output,
            contigs,
        } => vcf::run(
            bam,
            motifs,
            motif_file,
            sequence_fallback,
            sequence_index,
            output,
            sample_name,
            methylation_threshold,
            plain_output,
            contigs,
        ),
        Command::ImputeBam {
            bam,
            output,
            methylation_threshold,
            summary,
            motifs,
            motif_file,
            impute_all,
        } => impute_bam::run(
            bam,
            output,
            methylation_threshold,
            summary,
            motifs,
            motif_file,
            impute_all,
        ),
    }
}
