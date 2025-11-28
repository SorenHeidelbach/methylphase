use anyhow::Result;
use clap::Parser;
use methylphase::{cli::Cli, run};

fn main() -> Result<()> {
    let cli = Cli::parse();
    run(cli)
}
