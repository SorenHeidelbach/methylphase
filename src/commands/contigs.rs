use crate::io::BamReader;
use anyhow::{Context, Result};
use std::path::PathBuf;

pub fn run(bam: PathBuf) -> Result<()> {
    let reader =
        BamReader::open(&bam).with_context(|| format!("failed to open BAM {}", bam.display()))?;
    for contig in reader.contigs() {
        println!("{contig}");
    }
    Ok(())
}
