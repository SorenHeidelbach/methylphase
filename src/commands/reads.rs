use crate::io::BamReader;
use anyhow::{Context, Result};
use std::path::PathBuf;

pub fn run(bam: PathBuf, contigs: Vec<String>, limit: Option<usize>) -> Result<()> {
    let mut reader =
        BamReader::open(&bam).with_context(|| format!("failed to open BAM {}", bam.display()))?;
    let targets = if contigs.is_empty() {
        reader.contigs()
    } else {
        contigs
    };

    for target in targets {
        println!("# {target}");
        let reads = reader
            .read_ids(&target, limit)
            .with_context(|| format!("failed to query reads for {target}"))?;
        for read in reads {
            println!("{read}");
        }
    }

    Ok(())
}
