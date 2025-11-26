use crate::io::BamReader;
use anyhow::{Context, Result};
use std::{collections::HashSet, path::PathBuf};

pub fn run(bam: PathBuf, contigs: Vec<String>) -> Result<()> {
    let reader =
        BamReader::open(&bam).with_context(|| format!("failed to open BAM {}", bam.display()))?;
    let available = reader.contigs();
    if contigs.is_empty() {
        for contig in available {
            println!("{contig}");
        }
        return Ok(());
    }

    let available_set: HashSet<_> = available.into_iter().collect();
    for contig in contigs {
        if available_set.contains(&contig) {
            println!("{contig}");
        } else {
            eprintln!(
                "contigs: contig {} not present in {}; skipping",
                contig,
                bam.display()
            );
        }
    }

    Ok(())
}
