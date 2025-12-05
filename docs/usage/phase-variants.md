---
layout: default
title: Phase variants pipeline
---

# Phase variants pipeline

`methylphase phase-variants` is the end-to-end typing workflow. It extracts methylation features (via split-reads), merges them with a Floria haploset, selects the best latent class model, imputes labels, and emits summaries.

## Inputs
- `--floria <haploset>`: Floria haplosets to merge with methylation features. This file can be found in the floria output under `<out-dir>/<contig>/<contig>.haplosets`
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


## Full help
```
End-to-end methylation typing pipeline (formerly `methyltyping run`)

Usage: methylphase phase-variants [OPTIONS] --floria <FLORIA> --bam <BAM> --out <OUT>

Options:
  -h, --help  Print help

INPUT:
      --floria <FLORIA>           Floria haploset file
      --bam <BAM>                 BAM with methylation tags
      --sequence-fallback <FILE>  Optional FASTQ/BAM providing sequences when BAM entries omit SEQ (for split-reads)

MOTIFS:
      --motif [<SPEC>...]  Motif descriptors in format: motif_modtype_modposition. To specify multiple, separate with commas
      --motif-file <FILE>  Path to a motif TSV file. Should contain one motif descriptor per line with headers: motif, mod_type, mod_position

MODEL SELECTION:
      --min-classes <MIN_CLASSES>
          Minimum latent classes to evaluate (lower bound on number of phase variants) [default: 1]
      --max-classes <MAX_CLASSES>
          Maximum latent classes to evaluate (upper bound on number of phase variants) [default: 10]
      --criterion <CRITERION>
          Model selection criterion: bic (penalized likelihood), icl (bic + entropy), or cv (cross-val NLL) [default: icl] [possible values: bic, icl, cv]
      --penalty-multiplier <PENALTY_MULTIPLIER>
          Scale on the BIC/ICL penalty term; larger values favor fewer phase variants [default: 2.0]
      --cv-folds <CV_FOLDS>
          Folds for cross-validated scoring (used only when criterion=cv) [default: 5]

OUTPUT:
  -o, --out <OUT>  Output directory

PRIORS:
      --alpha-pi <ALPHA_PI>    Dirichlet prior on class proportions (higher -> smoother class weights) [default: 1.0]
      --alpha-phi <ALPHA_PHI>  Dirichlet prior on categorical emissions (higher -> smoother feature probabilities) [default: 1.0]

OPTIMIZATION:
      --max-iter <MAX_ITER>
          Maximum EM iterations per fit [default: 100]
      --tol <TOL>
          EM convergence tolerance on log-likelihood change [default: 1e-6]
      --min-class-weight <MIN_CLASS_WEIGHT>
          Drop classes whose average responsibility falls below this fraction; can roughly be interpreted as minimum phase variant abundance [default: 0.0]
```