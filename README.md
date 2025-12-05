# methylphase

Toolkit for BAM-driven methylation phasing and typing: extracts motif-level methylation from mod-tagged BAMs, clusters reads, merges with Floria haplotypes, and imputes labels with latent class models.

- Full [docs](https://sorenheidelbach.github.io/methylphase/) 

## Features
- End-to-end phasing (`phase-variants`) that joins methylation features with Floria haplosets and calls number of variants and assign reads.
- Split reads by methylation profiles; emit methylation cluster tables and optional FASTQ/BAM splits.
- Extract per-read and aggregate motif methylation tables.
- Utilities for exploring read level methylation profiles.

## Installation
Requires Rust (stable, 2021 edition) and a C toolchain.

Clone and install:
```bash
git clone https://github.com/SorenHeidelbach/methylphase.git
cd methylphase
cargo install --path . --force
methylphase --help
```

Or directly from GitHub:
```bash
cargo install --git https://github.com/SorenHeidelbach/methylphase --locked
```

## Quick start
From an indexed, mod-tagged BAM and a Floria haploset:

```bash
methylphase phase-variants \
  --floria haplotypes.hapset \
  --bam mapped_sample.mod.bam \
  --motif GATC_6mA_1 \
  --out phase_variants 
```

Key outputs: `split_reads/read_clustering_raw.tsv`, `split_reads/cluster_summary`, `dataset.tsv`, `imputed_labels.tsv`, `summary.json`.

Split-only and extraction workflow:

```bash
methylphase extract \
  --bam sample.mod.bam \
  --motif GATC_6mA_1 \
  --output-dir extract

methylphase split-reads \
  --bam sample.mod.bam \
  --motif GATC_6mA_1 \
  --output-dir split \
  --emit-fastq 
```

For arguments, inputs, and examples see the [GitHub Pages docs](https://sorenheidelbach.github.io/methylphase/).
