use crate::{core::MotifQuery, io::load_motif_file};
use anyhow::{bail, Result};
use std::path::PathBuf;

pub fn load_motif_queries(
    mut motifs: Vec<String>,
    motif_file: Option<PathBuf>,
) -> Result<Vec<MotifQuery>> {
    if let Some(path) = motif_file {
        let mut file_specs = load_motif_file(&path)?;
        motifs.append(&mut file_specs);
    }

    if motifs.is_empty() {
        bail!("no motifs provided; use --motif or --motif-file");
    }

    motifs.iter().map(|spec| MotifQuery::parse(spec)).collect()
}
