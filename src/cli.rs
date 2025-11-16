use clap::{Parser, Subcommand};
use std::path::PathBuf;

/// Command-line interface definition.
#[derive(Debug, Parser)]
#[command(
    name = "methylation_phasing",
    version,
    about = "Tools for BAM-driven methylation phasing"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// List contigs present in the BAM.
    Contigs {
        /// Path to the indexed BAM file.
        #[arg(value_name = "BAM")]
        bam: PathBuf,
    },
    /// Emit read identifiers for the requested contigs.
    Reads {
        /// Path to the indexed BAM file.
        #[arg(value_name = "BAM")]
        bam: PathBuf,

        /// Contigs to query. If omitted all contigs are used.
        #[arg(short, long, value_name = "CONTIG", num_args = 1.., value_delimiter = ',')]
        contigs: Vec<String>,

        /// Maximum number of read identifiers to emit per contig.
        #[arg(long)]
        limit: Option<usize>,
    },
    /// Extract methylation calls that overlap the requested motifs.
    Extract {
        /// Path to the indexed BAM file.
        #[arg(value_name = "BAM")]
        bam: PathBuf,

        /// Motif descriptors in the format motif_modtype_modposition.
        #[arg(
            short = 'm',
            long = "motif",
            value_name = "SPEC",
            num_args = 1..,
            value_delimiter = ','
        )]
        motifs: Vec<String>,

        /// Path to a TSV/CSV with columns motif, mod_type, mod_position.
        #[arg(long = "motif-file", value_name = "FILE")]
        motif_file: Option<PathBuf>,

        /// Optional TSV where each row reports if a read is methylated for each motif.
        #[arg(long = "motif-summary-tsv", value_name = "FILE")]
        motif_summary_tsv: Option<PathBuf>,

        /// Optional directory to write FASTQ subsets grouped by methylated motif combinations.
        #[arg(long = "fastq-dir", value_name = "DIR")]
        fastq_dir: Option<PathBuf>,

        /// Path to the per-read TSV output.
        #[arg(long = "per-read-tsv", value_name = "FILE")]
        per_read_tsv: PathBuf,

        /// Path to the per-read aggregate TSV output.
        #[arg(long = "aggregate-tsv", value_name = "FILE")]
        aggregate_tsv: PathBuf,

        /// Restrict processing to the provided contigs.
        #[arg(
            short = 'c',
            long = "contig",
            value_name = "CONTIG",
            value_delimiter = ','
        )]
        contigs: Vec<String>,
    },
    /// Evaluate longitudinal methylation agreement across contigs.
    Longitudinal {
        /// Path to the indexed BAM file.
        #[arg(value_name = "BAM")]
        bam: PathBuf,

        /// Motif descriptors in the format motif_modtype_modposition.
        #[arg(
            short = 'm',
            long = "motif",
            value_name = "SPEC",
            num_args = 1..,
            value_delimiter = ','
        )]
        motifs: Vec<String>,

        /// Path to a TSV/CSV with columns motif, mod_type, mod_position.
        #[arg(long = "motif-file", value_name = "FILE")]
        motif_file: Option<PathBuf>,

        /// Path to the TSV output describing per-block read pattern frequencies.
        #[arg(long = "output", value_name = "FILE")]
        output: PathBuf,

        /// Block size in base pairs used to aggregate motif positions.
        #[arg(long = "block-size", value_name = "BP", default_value_t = 1000)]
        block_size: usize,

        /// Minimum probability to classify a modification call as methylated.
        #[arg(long = "methylation-threshold", default_value_t = 0.5)]
        methylation_threshold: f32,

        /// Restrict processing to the provided contigs.
        #[arg(
            short = 'c',
            long = "contig",
            value_name = "CONTIG",
            value_delimiter = ','
        )]
        contigs: Vec<String>,
    },
}
