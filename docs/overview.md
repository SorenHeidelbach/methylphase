# Overview

methylphase is a toolkit for BAM-driven methylation phasing and typing. It extracts per-read methylation features from mod-tagged alignments (MM/ML tags), clusters reads by motif patterns, and aggregates or imputes results for downstream analyses.

## What methylphase does
- **Phase variants pipeline**: end-to-end typing that merges motif-derived methylation features with a Floria haploset, selects a latent class model, and emits imputed labels and summaries.
- **Split reads**: cluster reads by motif methylation profiles (GMM by default; HDBSCAN and agglomerative are available) and write clustering tables plus optional FASTQ/BAM splits.
- **Extract**: intersect motifs with alignments to produce per-read and aggregate methylation tables, and optionally motif summaries or FASTQ subsets.
- **Typing utilities**: EM fitting, model selection, imputation, Floria conversion, and FASTQ splitting from existing assignments.
- **Utilities**: contig/read inspection, VCF summarization, methylation imputation directly into BAM, and longitudinal agreement checks.

## Inputs and assumptions
- **Aligned BAM with MM/ML tags** is required for methylation coordinates; provide a BAM index. If SEQ is missing, supply `--sequence-fallback` (FASTQ/BAM).
- **Motifs** can be inline (`GATC_6mA_1`) or provided via file (TSV/CSV with motif, mod_type, mod_position, optional complement/bin/reference columns).
- **Floria haploset** (for `phase-variants`) supplies SNP haplotypes that are merged with methylation features.
- Motif interpretation is strand-aware; complements in motif files control reverse-strand handling.

## Outputs (high level)
- **Split reads**: `read_clustering.tsv`, `read_clustering_raw.tsv`, `cluster_summary.tsv`, and optional clustered FASTQ/BAMs.
- **Phase variants**: internal `split_reads/` artifacts, `dataset.tsv`, `categories.toml`, `fits/`, `best_model.json`, `best_responsibilities.tsv`, `imputed.tsv`, `imputed_labels.tsv`, and `summary.json`.
- **Extract**: `per_read.tsv`, `aggregate.tsv`, optional `motif_summary.tsv`, and FASTQ subsets when requested.
- **Utilities**: VCFs (`*.vcf` or bgzip+tabix), imputed BAMs, and longitudinal reports.

## Modeling notes
- Latent class model uses EM with Dirichlet priors for weights/emissions; supports continuous methylation features (Gaussian) and categorical haplotypes.
- Model selection defaults to ICL with penalty multiplier 2.0 and class bounds 1..10; alternative criteria include BIC and cross-validation.
- Imputation outputs MAP labels for categorical fields and responsibility-weighted means for methylation features; classes with low weight can be pruned.

## Project layout (code)
- `src/main.rs`: CLI entrypoint, delegates to `commands::run`.
- `src/cli.rs`: top-level CLI definitions for commands/groups.
- `src/commands/`: BAM extraction, split-reads clustering, VCF export, BAM imputation, contig/read utilities, longitudinal agreement, and shared helpers.
- `src/typing/`: typing entrypoint/CLI plus data handling, EM, model selection, imputation, methylation feature logic, Floria parsing, FASTQ splitting, and priors.

See the usage pages for command-by-command details and examples.
