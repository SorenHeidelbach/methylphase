use anyhow::Result;
use clap::Parser;
use methylation_phasing::{cli::Cli, run};

fn main() -> Result<()> {
    let cli = Cli::parse();
    run(cli)
}
