pub mod contigs;
pub mod extract;
pub mod longitudinal;
pub mod reads;
pub mod shared;

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
            per_read_tsv,
            aggregate_tsv,
            contigs,
        } => extract::run(
            bam,
            motifs,
            motif_file,
            motif_summary_tsv,
            fastq_dir,
            per_read_tsv,
            aggregate_tsv,
            contigs,
        ),
        Command::Longitudinal {
            bam,
            motifs,
            motif_file,
            output,
            block_size,
            methylation_threshold,
            contigs,
        } => longitudinal::run(
            bam,
            motifs,
            motif_file,
            output,
            block_size,
            methylation_threshold,
            contigs,
        ),
    }
}
