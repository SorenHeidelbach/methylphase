---
layout: default
title: Quickstart
---

# Quickstart

Use this minimal workflow to go from an indexed, mod-tagged BAM to phased methylation labels.

## 1) Collect inputs
- Indexed BAM with MM/ML tags: `sample.mod.bam` and `sample.mod.bam.bai`.
- Motifs to interrogate (inline list or motif file).
- Floria haploset (for the full typing pipeline).

## 2) Run the end-to-end pipeline
This drives split-reads internally, merges haplotype and methylation features, selects the best latent class model, and imputes labels.
```bash
methylphase phase-variants \
  --floria haplotypes.hapset \
  --bam sample.mod.bam \
  --motif CG_5,GATC_6mA_1 \
  --out results/phase_variants \
  --min-classes 1 \
  --max-classes 8
```
Key outputs land under `results/phase_variants`: `split_reads/*`, `dataset.tsv`, `categories.toml`, `best_model.json`, `best_responsibilities.tsv`, `imputed.tsv`, and `imputed_labels.tsv`.

## 3) Inspect quick summaries
- Cluster assignments and raw features: `split_reads/read_clustering_raw.tsv`
- Imputed labels and summary JSON: `imputed_labels.tsv`, `summary.json`
- Model diagnostics per class count: `fits/`

## Alternative quick run: split and aggregate only
If you just need read clustering and per-read methylation tables:
```bash
methylphase extract \
  --bam sample.mod.bam \
  --motif CG_5 \
  --output-dir results/extract

methylphase split-reads \
  --bam sample.mod.bam \
  --motif-file motifs.tsv \
  --output-dir results/split \
  --cluster-algorithm hdbscan \
  --emit-fastq \
  --threads 8
```

Next, dive into the command-specific pages under `usage/` for detailed arguments and outputs.
