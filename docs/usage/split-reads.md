---
layout: default
title: Split reads
---

# Split reads

`methylphase split-reads` clusters reads by motif methylation features and can write clustering tables plus per-cluster FASTQ/BAM files.

## Common arguments
- `--bam <mod.bam>`: indexed BAM with MM/ML tags.
- Motifs: `--motif ...` or `--motif-file ...` (same format as `extract`).
- `--output-dir <dir>`: destination for clustering outputs.
- `--cluster-algorithm <gmm|hdbscan|agglomerative>`: clustering backend (default `gmm`).
- `--threads <n>`: parallelism for feature extraction and clustering.
- Optional `--emit-fastq` (and related flags) to write per-cluster FASTQ/BAMs.

## Outputs
- `read_clustering_raw.tsv`: per-read cluster assignments with raw motif features.
- `read_clustering.tsv`: cleaned/filtered clustering table.
- `cluster_summary.tsv`: per-cluster counts and feature summaries.
- Optional `clusters/*.fastq` or `*.bam` when emitting sequences.

## Example
```bash
methylphase split-reads \
  --bam sample.mod.bam \
  --motif CG_5,GATC_6mA_1 \
  --output-dir results/split \
  --cluster-algorithm hdbscan \
  --emit-fastq \
  --threads 8
```
