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

        /// Motif descriptors in the format motif_modtype_modposition (0-based modification index).
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

        /// Optional FASTQ/BAM providing sequences for reads whose BAM entries omit SEQ.
        #[arg(long = "sequence-fallback", value_name = "FILE")]
        sequence_fallback: Option<PathBuf>,

        /// Optional index (e.g., fqidx) enabling random access to the sequence fallback FASTQ.
        #[arg(long = "sequence-index", value_name = "FILE")]
        sequence_index: Option<PathBuf>,

        /// Directory to place all outputs (per-read.tsv, aggregate.tsv, etc.).
        #[arg(long = "output-dir", value_name = "DIR")]
        output_dir: PathBuf,

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

        /// Motif descriptors in the format motif_modtype_modposition (0-based modification index).
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

        /// Directory to place the longitudinal output (longitudinal.tsv).
        #[arg(long = "output-dir", value_name = "DIR")]
        output_dir: PathBuf,

        /// Block size in base pairs used to aggregate motif positions.
        #[arg(long = "block-size", value_name = "BP", default_value_t = 1000)]
        block_size: usize,

        /// Minimum probability to classify a modification call as methylated.
        #[arg(long = "methylation-threshold", default_value_t = 0.5)]
        methylation_threshold: f32,

        /// Optional FASTQ/BAM providing sequences for reads whose BAM entries omit SEQ.
        #[arg(long = "sequence-fallback", value_name = "FILE")]
        sequence_fallback: Option<PathBuf>,

        /// Optional index (e.g., fqidx) enabling random access to the sequence fallback FASTQ.
        #[arg(long = "sequence-index", value_name = "FILE")]
        sequence_index: Option<PathBuf>,

        /// Restrict processing to the provided contigs.
        #[arg(
            short = 'c',
            long = "contig",
            value_name = "CONTIG",
            value_delimiter = ','
        )]
        contigs: Vec<String>,
    },
    /// Cluster reads by motif methylation and split FASTQs by cluster.
    SplitReads {
        /// Path to the indexed BAM file.
        #[arg(value_name = "BAM")]
        bam: PathBuf,

        /// Motif descriptors in the format motif_modtype_modposition (0-based modification index).
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

        /// Optional FASTQ/BAM providing sequences for reads whose BAM entries omit SEQ.
        #[arg(long = "sequence-fallback", value_name = "FILE")]
        sequence_fallback: Option<PathBuf>,

        /// Optional index (e.g., fqidx) enabling random access to the sequence fallback FASTQ.
        #[arg(long = "sequence-index", value_name = "FILE")]
        sequence_index: Option<PathBuf>,

        /// Directory to place outputs (per_read.tsv, aggregate.tsv, read_clustering.tsv, clusters/).
        #[arg(long = "output-dir", value_name = "DIR")]
        output_dir: PathBuf,

        /// Minimum cluster size for HDBSCAN (acts as a floor; the tool auto-scales upward when many reads are available).
        #[arg(long = "min-cluster-size", default_value_t = 5)]
        min_cluster_size: usize,

        /// Minimum samples for density estimation (defaults to min-cluster-size if omitted).
        #[arg(long = "min-samples")]
        min_samples: Option<usize>,

        /// Emit FASTQ cluster outputs instead of BAM.
        #[arg(long = "emit-fastq")]
        emit_fastq: bool,

        /// Number of worker threads to use when processing reads.
        #[arg(long = "threads", default_value_t = 1)]
        threads: usize,

        /// Restrict processing to the provided contigs.
        #[arg(
            short = 'c',
            long = "contig",
            value_name = "CONTIG",
            value_delimiter = ','
        )]
        contigs: Vec<String>,
    },
    /// Summarize motif methylation as a VCF.
    Vcf {
        /// Path to the indexed BAM file.
        #[arg(value_name = "BAM")]
        bam: PathBuf,

        /// Motif descriptors in the format motif_modtype_modposition (0-based modification index).
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

        /// Optional FASTQ/BAM providing sequences for reads whose BAM entries omit SEQ.
        #[arg(long = "sequence-fallback", value_name = "FILE")]
        sequence_fallback: Option<PathBuf>,

        /// Optional index (e.g., fqidx) enabling random access to the sequence fallback FASTQ.
        #[arg(long = "sequence-index", value_name = "FILE")]
        sequence_index: Option<PathBuf>,

        /// Output VCF path.
        #[arg(long = "output", value_name = "VCF")]
        output: PathBuf,

        /// Sample name to write into the VCF (defaults to \"sample\").
        #[arg(long = "sample-name", value_name = "NAME")]
        sample_name: Option<String>,

        /// Minimum probability to classify a modification call as methylated.
        #[arg(long = "methylation-threshold", default_value_t = 0.5)]
        methylation_threshold: f32,

        /// Disable bgzip/tabix generation and leave the VCF uncompressed.
        #[arg(long = "plain-output")]
        plain_output: bool,

        /// Restrict processing to the provided contigs.
        #[arg(
            short = 'c',
            long = "contig",
            value_name = "CONTIG",
            value_delimiter = ','
        )]
        contigs: Vec<String>,
    },
    /// Impute high-confidence methylation calls directly into the BAM.
    ImputeBam {
        /// Path to the source BAM file (must include MM/ML tags).
        #[arg(value_name = "BAM")]
        bam: PathBuf,

        /// Path to the output BAM file with imputed bases.
        #[arg(long = "output", value_name = "BAM")]
        output: PathBuf,

        /// Minimum probability to imprint a modification as a base/quality change.
        #[arg(long = "methylation-threshold", default_value_t = 0.8)]
        methylation_threshold: f32,

        /// Optional path for the imputation summary TSV (defaults to output with .summary.tsv).
        #[arg(long = "summary", value_name = "FILE")]
        summary: Option<PathBuf>,

        /// Motif descriptors in the format motif_modtype_modposition (0-based modification index).
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

        /// Impute every high-confidence modification (ignores motifs).
        #[arg(long = "impute-all", conflicts_with_all = ["motifs", "motif_file"])]
        impute_all: bool,
    },
}
