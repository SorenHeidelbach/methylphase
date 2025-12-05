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


## Full help
```
Cluster reads by motif methylation and split FASTQs by cluster

Usage: methylphase split-reads [OPTIONS] --output-dir <DIR> <BAM>

Arguments:
  <BAM>  Path to the indexed BAM file

Options:
  -m, --motif <SPEC>...
          Motif descriptors in the format motif_modtype_modposition (0-based modification index)
      --motif-file <FILE>
          Path to a TSV/CSV with columns motif, mod_type, mod_position (optional motif_complement, mod_position_complement, and id/bin/reference for bin-specific motifs)
      --sequence-fallback <FILE>
          Optional FASTQ/BAM providing sequences for reads whose BAM entries omit SEQ (indexes stored as <fallback>.fqidx/.bai)
      --sequence-index <FILE>
          Optional override for the fallback index path (defaults to <fallback>.fqidx/.bai)
      --output-dir <DIR>
          Directory to place outputs (per_read.tsv, aggregate.tsv, read_clustering.tsv, clusters/)
      --min-cluster-size <MIN_CLUSTER_SIZE>
          Minimum cluster size for clustering (acts as a floor; auto-scales upward when many reads are available) [default: 5]
      --min-samples <MIN_SAMPLES>
          Minimum samples for density estimation (defaults to min-cluster-size if omitted)
      --emit-fastq
          Emit FASTQ cluster outputs instead of BAM
      --threads <THREADS>
          Number of worker threads to use when processing reads [default: 1]
      --cluster-algorithm <CLUSTER_ALGORITHM>
          Clustering algorithm to use for assigning reads to groups [default: gmm] [possible values: hdbscan, gmm, agglomerative]
  -c, --contig <CONTIG>
          Restrict processing to the provided contigs
      --contig-bins <TSV>
          Optional TSV/CSV mapping contig identifiers to bins (contig<TAB>bin)
  -h, --help
          Print help
```