use thiserror::Error;

/// Custom error type for the CLI and library.
#[derive(Debug, Error)]
pub enum MethylError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error: {0}")]
    Parse(String),
    #[error("Config error: {0}")]
    Config(String),
    #[error("CSV error: {0}")]
    Csv(#[from] csv::Error),
    #[error("TOML error: {0}")]
    Toml(#[from] toml::de::Error),
    #[error("Inconsistent category {category}: {message}")]
    InconsistentCategory { category: usize, message: String },
    #[error("Shape error: {0}")]
    Shape(String),
    #[error("Math error: {0}")]
    Math(String),
    #[error("Other: {0}")]
    Other(String),
}
