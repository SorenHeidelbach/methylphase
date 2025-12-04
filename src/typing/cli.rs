use clap::{ArgGroup, Subcommand};
use std::path::PathBuf;

/// Supported subcommands.
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// End-to-end pipeline: parse Floria + methylation, select best class count, impute labels.
    #[command(
        about = "End-to-end methylation typing pipeline (runs split-reads internally if needed)",
        long_about = "Runs the full methylation typing pipeline: combines a Floria haploset with \
read-level methylation features (either precomputed or generated via split-reads), selects the \
best latent class count, imputes missing values, and emits labels and artifacts under --out."
    )]
    Run(RunArgs),
    /// Fit an EM model with a fixed number of classes.
    Fit(FitArgs),
    /// Select best number of classes using BIC.
    SelectBestC(SelectArgs),
    /// Impute missing entries using a fitted model.
    Impute(ImputeArgs),
    /// Convert a Floria haploset file into dataset/config TSV+TOML.
    ConvertFloria(ConvertFloriaArgs),
    /// Split FASTQ into per-class FASTQ files using assignments.
    SplitFastq(SplitFastqArgs),
}

/// Arguments for fitting a model.
#[derive(clap::Args, Debug)]
pub struct FitArgs {
    #[arg(long)]
    pub data: PathBuf,
    #[arg(long)]
    pub config: PathBuf,
    #[arg(long)]
    pub classes: usize,
    #[arg(long, default_value = "1.0")]
    pub alpha_pi: f64,
    #[arg(long, default_value = "1.0")]
    pub alpha_phi: f64,
    #[arg(long, default_value = "100")]
    pub max_iter: usize,
    #[arg(long, default_value = "1e-6")]
    pub tol: f64,
    /// Minimum effective class weight (fraction of data) to keep after EM.
    #[arg(long, default_value = "0.0")]
    pub min_class_weight: f64,
    /// Delimiter for input data. Default is tab.
    #[arg(long, default_value = "\t")]
    pub delimiter: String,
    #[arg(long)]
    pub output_model: PathBuf,
    /// Optional output for responsibilities (per-block posterior).
    #[arg(long)]
    pub output_responsibilities: Option<PathBuf>,
}

/// Arguments for selecting number of classes.
#[derive(clap::Args, Debug)]
pub struct SelectArgs {
    #[arg(long)]
    pub data: PathBuf,
    #[arg(long)]
    pub config: PathBuf,
    #[arg(long)]
    pub min_classes: usize,
    #[arg(long)]
    pub max_classes: usize,
    #[arg(long, default_value = "1.0")]
    pub alpha_pi: f64,
    #[arg(long, default_value = "1.0")]
    pub alpha_phi: f64,
    /// Model selection criterion.
    #[arg(long, value_parser = ["bic", "icl", "cv"], default_value = "icl")]
    pub criterion: String,
    /// Multiplier applied to the BIC/ICL penalty term to discourage extra classes.
    #[arg(long, default_value = "2.0")]
    pub penalty_multiplier: f64,
    /// Number of folds for cross-validated selection (criterion=cv).
    #[arg(long, default_value = "5")]
    pub cv_folds: usize,
    #[arg(long, default_value = "100")]
    pub max_iter: usize,
    #[arg(long, default_value = "1e-6")]
    pub tol: f64,
    /// Minimum effective class weight (fraction of data) to keep after EM.
    #[arg(long, default_value = "0.0")]
    pub min_class_weight: f64,
    #[arg(long, default_value = "\t")]
    pub delimiter: String,
    #[arg(long)]
    pub output_model: PathBuf,
}

/// Arguments for imputation.
#[derive(clap::Args, Debug)]
pub struct ImputeArgs {
    #[arg(long)]
    pub data: PathBuf,
    #[arg(long)]
    pub config: PathBuf,
    #[arg(long)]
    pub model: PathBuf,
    #[arg(long)]
    pub output: PathBuf,
    #[arg(long, default_value = "\t")]
    pub delimiter: String,
}

/// Arguments for converting Floria haploset output.
#[derive(clap::Args, Debug)]
pub struct ConvertFloriaArgs {
    /// Path to Floria haploset file.
    #[arg(long)]
    pub floria: PathBuf,
    /// Optional methylation cluster assignments TSV (read_id, cluster_id,...).
    #[arg(long)]
    pub methylation: Option<PathBuf>,
    /// Output dataset TSV/CSV.
    #[arg(long)]
    pub output_data: PathBuf,
    /// Output category config TOML.
    #[arg(long)]
    pub output_config: PathBuf,
    /// Delimiter for output data.
    #[arg(long, default_value = "\t")]
    pub delimiter: String,
}

/// Arguments for the full pipeline.
#[derive(clap::Args, Debug)]
pub struct RunArgs {
    /// Floria haploset file.
    #[arg(long, help_heading = "INPUT")]
    pub floria: PathBuf,
    /// BAM with methylation tags.
    #[arg(long, help_heading = "INPUT")]
    pub bam: PathBuf,
    /// Motif descriptors in format: motif_modtype_modposition, e.g. GATC_6mA_1. To specify multiple separate with commas, e.g. GATC_6mA_1,CCWGG_5mC_1.
    #[arg(
        long = "motif",
        value_name = "SPEC",    
        num_args = 0..,
        value_delimiter = ',',
        help_heading = "MOTIFS"
    )]
    pub motifs: Vec<String>,
    /// Path to a motif TSV file. Should contain one motif descriptor per line with headers: motif, mod_type, mod_position.
    #[arg(long = "motif-file", value_name = "FILE", help_heading = "MOTIFS")]
    pub motif_file: Option<PathBuf>,
    /// Optional FASTQ/BAM providing sequences when BAM entries omit SEQ.
    #[arg(long = "sequence-fallback", value_name = "FILE", help_heading = "INPUT")]
    pub sequence_fallback: Option<PathBuf>,
    // split-reads tuning uses defaults; use the split-reads command for manual control.
    /// Minimum latent classes to evaluate (lower bound on number of phase variants).
    #[arg(long, default_value_t = 1, help_heading = "MODEL SELECTION")]
    pub min_classes: usize,
    /// Maximum latent classes to evaluate (upper bound on number of phase variants).
    #[arg(long, default_value_t = 10, help_heading = "MODEL SELECTION")]
    pub max_classes: usize,
    /// Output directory.
    #[arg(short = 'o', long, help_heading = "OUTPUT")]
    pub out: PathBuf,
    /// Dirichlet prior on class proportions (higher -> smoother class weights).
    #[arg(long, default_value = "1.0", help_heading = "PRIORS")]
    pub alpha_pi: f64,
    /// Dirichlet prior on categorical emissions (higher -> smoother feature probabilities).
    #[arg(long, default_value = "1.0", help_heading = "PRIORS")]
    pub alpha_phi: f64,
    /// Model selection criterion: bic (penalized likelihood), icl (bic + entropy), or cv (cross-val NLL).
    #[arg(
        long,
        value_parser = ["bic", "icl", "cv"],
        default_value = "icl",
        help_heading = "MODEL SELECTION"
    )]
    pub criterion: String,
    /// Scale on the BIC/ICL penalty term; larger values favor fewer phase variants.
    #[arg(
        long,
        default_value = "2.0",
        help_heading = "MODEL SELECTION"
    )]
    pub penalty_multiplier: f64,
    /// Folds for cross-validated scoring (used only when criterion=cv).
    #[arg(
        long,
        default_value = "5",
        help_heading = "MODEL SELECTION"
    )]
    pub cv_folds: usize,
    /// Maximum EM iterations per fit.
    #[arg(long, default_value = "100", help_heading = "OPTIMIZATION")]
    pub max_iter: usize,
    /// EM convergence tolerance on log-likelihood change.
    #[arg(long, default_value = "1e-6", help_heading = "OPTIMIZATION")]
    pub tol: f64,
    /// Drop classes whose average responsibility falls below this fraction; can roughly be interpreted as minimum phase variant abundance.
    #[arg(
        long,
        default_value = "0.0",
        help_heading = "OPTIMIZATION"
    )]
    pub min_class_weight: f64,
    #[arg(long, default_value = "\t", help_heading = "IO")]
    pub delimiter: String,
}

/// Arguments for splitting FASTQ by class assignments.
#[derive(clap::Args, Debug)]
pub struct SplitFastqArgs {
    /// FASTQ file to split.
    #[arg(long)]
    pub fastq: PathBuf,
    /// Assignments TSV (imputed_labels or responsibilities).
    #[arg(long)]
    pub assignments: PathBuf,
    /// Output directory for per-class FASTQ files.
    #[arg(short = 'o', long)]
    pub out: PathBuf,
    #[arg(long, default_value = "\t")]
    pub delimiter: String,
}

impl FitArgs {
    /// Resolve delimiter to a single byte.
    pub fn delimiter_byte(&self) -> u8 {
        self.delimiter.as_bytes()[0]
    }
}

impl SelectArgs {
    pub fn delimiter_byte(&self) -> u8 {
        self.delimiter.as_bytes()[0]
    }
}

impl ImputeArgs {
    pub fn delimiter_byte(&self) -> u8 {
        self.delimiter.as_bytes()[0]
    }
}

impl ConvertFloriaArgs {
    pub fn delimiter_byte(&self) -> u8 {
        self.delimiter.as_bytes()[0]
    }
}

impl RunArgs {
    pub fn delimiter_byte(&self) -> u8 {
        self.delimiter.as_bytes()[0]
    }
}

impl SplitFastqArgs {
    pub fn delimiter_byte(&self) -> u8 {
        self.delimiter.as_bytes()[0]
    }
}
