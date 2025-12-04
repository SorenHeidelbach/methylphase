# Installation

## Requirements
- Rust (stable, 2021 edition) and a C toolchain for compiling Rust crates.
- Indexed BAM inputs with MM/ML tags and available disk space for outputs.
- Optional: `bgzip` and `tabix` if you plan to bgzip/tabix VCFs produced by `utils vcf`.

## Install from source
```bash
git clone <this-repo>
cd methylphase
cargo build --release           # produces target/release/methylphase
# or install into ~/.cargo/bin:
cargo install --path . --force
methylphase --help              # verify the binary runs
```

If you prefer a one-liner without cloning, you can install directly from the repository URL once it is public:  
`cargo install --git https://github.com/<org>/methylphase --locked`

## Updating
Re-run `cargo install --path . --force` (or the `--git` form) after pulling changes.

## Environment tips
- Ensure BAM files are indexed (`samtools index sample.bam`) before running commands.
- Provide `--sequence-fallback` when SEQ is absent in BAM headers (FASTQ/BAM supported).
- Motifs can be inline (`GATC_6mA_1`) or supplied via TSV/CSV; see [extract](../usage/extract.md) for details.
