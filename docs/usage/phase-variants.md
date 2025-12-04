# Phase variants pipeline

`methylphase phase-variants` is the end-to-end typing workflow. It extracts methylation features (via split-reads), merges them with a Floria haploset, selects the best latent class model, imputes labels, and emits summaries.

## Inputs
- `--floria <haploset>`: Floria haplotype set to merge with methylation features.
- `--bam <mod.bam>`: Indexed BAM with MM/ML tags.
- Motifs: `--motif CG_5,GATC_6mA_1` or `--motif-file motifs.tsv`.
- Optional `--sequence-fallback <fastq|bam>` when SEQ is missing from the BAM.
- Optional binning (`--binning-file`) to stratify outputs by contig bin.

## Core options
- `--min-classes`, `--max-classes`: bounds for class search (default 1..10).
- `--criterion`: model selection (`icl` default, `bic`, or `cv`).
- `--penalty-multiplier`: scales the selection penalty (default 2.0 for ICL).
- `--out`: output directory (created if missing).
- Split-reads internals: uses GMM clustering by default to generate `read_clustering_raw.tsv`.

## Typical command
```bash
methylphase phase-variants \
  --floria haplotypes.hapset \
  --bam sample.mod.bam \
  --motif CG_5,GATC_6mA_1 \
  --out results/phase_variants \
  --min-classes 1 \
  --max-classes 8
```

## Outputs
- `split_reads/`: internal clustering artifacts, including `read_clustering_raw.tsv`.
- `dataset.tsv`, `categories.toml`: merged data/config used for typing.
- `fits/`: model fits per class count.
- `best_model.json`, `best_responsibilities.tsv`: selected model and responsibilities.
- `imputed.tsv`, `imputed_labels.tsv`: completed data and labels.
- `summary.json`: run metadata and selection summary.
