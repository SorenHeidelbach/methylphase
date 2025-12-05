---
layout: default
title: Installation
---

# Installation

## Requirements
- Rust (stable, 2021 edition) and a C toolchain for compiling Rust crates.
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

If you prefer a one-liner without cloning, you can install directly from the repository URL:  
`cargo install --git https://github.com/SorenHeidelbach/methylphase --locked`

