use crate::{core::MotifQuery, io::load_motif_file};
use anyhow::{bail, Result};
use std::path::PathBuf;

pub fn load_motif_queries(
    mut motifs: Vec<String>,
    motif_file: Option<PathBuf>,
    bin_ids: Option<&[String]>,
) -> Result<Vec<MotifQuery>> {
    if let Some(path) = motif_file {
        let mut file_specs = load_motif_file(&path, bin_ids)?;
        motifs.append(&mut file_specs);
    }

    if motifs.is_empty() {
        bail!("no motifs provided; use --motif or --motif-file");
    }

    eprintln!(
        "methylphase: loaded {} motif(s): {}",
        motifs.len(),
        motifs.join(",")
    );

    motifs.iter().map(|spec| MotifQuery::parse(spec)).collect()
}
