---
layout: default
title: Typing utilities
---

# Typing utilities

`methylphase typing <subcommand>` provides subcommands used in phase-variants. These can be used to fit models, impute missing values, convert Floria haplosets, and split FASTQ files by assigned classes. These command are mostly for internal use but may be helpful for users.

## fit
Fit a latent class model with a fixed number of classes.
```bash
methylphase typing fit \
  --data dataset.tsv \
  --config categories.toml \
  --classes 4 \
  --alpha-pi 1.0 \
  --alpha-phi 1.0 \
  --output-model fits/model_c4.json
```
Outputs a model JSON and (optionally) responsibilities.

## select-best-c
Score multiple class counts and pick the best according to a criterion.
```bash
methylphase typing select-best-c \
  --data dataset.tsv \
  --config categories.toml \
  --classes 1..10 \
  --criterion icl \
  --penalty-multiplier 2.0 \
  --output-dir fits
```
Emits per-class fits in `fits/`, `best_model.json`, and `best_responsibilities.tsv`.

## impute
Impute missing categorical/methylation values given a fitted model.
```bash
methylphase typing impute \
  --data dataset.tsv \
  --config categories.toml \
  --model best_model.json \
  --responsibilities best_responsibilities.tsv \
  --output imputed.tsv \
  --output-labels imputed_labels.tsv
```

## convert-floria
Convert a Floria haploset to the dataset/config format used by typing.
```bash
methylphase typing convert-floria \
  --haploset haplotypes.hapset \
  --output-data dataset.tsv \
  --output-config categories.toml
```
Optionally merge methylation features if provided.

## split-fastq
Split FASTQ files using assignments/responsibilities.
```bash
methylphase typing split-fastq \
  --fastq reads.fastq.gz \
  --responsibilities best_responsibilities.tsv \
  --output-dir fastq_by_class
```
