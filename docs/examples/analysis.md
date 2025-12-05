---
layout: default
title: Example analysis
---

# Example analysis

This walkthrough demonstrates a typical methylphase run on a mod-tagged BAM with motif-driven typing.

## Inputs
- `sample.mod.bam` + index
- Motifs file `motifs.tsv` (columns: `motif`, `mod_type`, `mod_position`, optional `complement`, `bin`, `reference`)
- Floria haploset `haplotypes.hapset`

## Steps
1) **(Optional) Inspect contigs**
   ```bash
   methylphase utils contigs --bam sample.mod.bam
   ```
2) **(Optional) Extract per-read methylation**
   ```bash
   methylphase extract \
     --bam sample.mod.bam \
     --motif-file motifs.tsv \
     --output-dir extract \
     --motif-summary
   ```
3) **Run the full phase-variants pipeline**
   ```bash
   methylphase phase-variants \
     --floria floria/contig/contig.haplosets \
     --bam sample.mod.bam \
     --motif-file motifs.tsv \
     --out phase_variants 
   ```
   
4) **Review clustering and model choice**

   - `phase_variants/split_reads/read_clustering_raw.tsv` for per-read motif features and clusters.
   - `phase_variants/fits/` for per-class log-likelihoods and criteria.
   - `phase_variants/summary.json` for the selected model and criteria values.

5) **Use the imputed outputs**

   - `imputed.tsv` for completed datasets (continuous and categorical fields).
   - `imputed_labels.tsv` for MAP class labels.
   - Feed these into downstream association tests or visualization tools.

6) **Split FASTQs by inferred class for downstream assembly/analysis of each variant:**
    ```bash
    methylphase typing split-fastq \
      --fastq reads.fastq \
      --responsibilities phase_variants/best_responsibilities.tsv \
      --output-dir phase_variants/fastq_by_class
    ```
