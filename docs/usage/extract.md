---
layout: default
title: Extract methylation calls
---

# Extract methylation calls

`methylphase extract` intersects motifs with mod-tagged alignments and emits per-read and aggregate methylation calls.

## Common arguments
- `--bam <mod.bam>`: indexed BAM with MM/ML tags.
- Motifs: `--motif CG_5mC_0,GATC_6mA_1` or `--motif-file motifs.tsv` (TSV/CSV with `motif`, `mod_type`, `mod_position`).
- `--output-dir <dir>`: destination for TSVs (created if missing).
- `--sequence-fallback <fastq|bam>`: when SEQ is absent in BAM.


## Outputs
- `per_read.tsv`: methylation calls per read and motif occurrence.
- `aggregate.tsv`: aggregated motif-level calls.

## Example
```bash
methylphase extract \
  --bam sample.mod.bam \
  --motif-file motifs.tsv \
  --output-dir results/extract 
```

## Full help
```
Extract methylation calls that overlap the requested motifs

Usage: methylphase extract [OPTIONS] --output-dir <DIR> <BAM>

Arguments:
  <BAM>  Path to the indexed BAM file

Options:
  -m, --motif <SPEC>...           Motif descriptors in the format motif_modtype_modposition (0-based modification index)
      --motif-file <FILE>         Path to a TSV/CSV with columns motif, mod_type, mod_position (optional motif_complement, mod_position_complement, and id/bin/reference for bin-specific motifs)
      --motif-summary-tsv <FILE>  Optional TSV where each row reports if a read is methylated for each motif
      --fastq-dir <DIR>           Optional directory to write FASTQ subsets grouped by methylated motif combinations
      --sequence-fallback <FILE>  Optional FASTQ/BAM providing sequences for reads whose BAM entries omit SEQ (indexes stored as <fallback>.fqidx/.bai)
      --sequence-index <FILE>     Optional override for the fallback index path (defaults to <fallback>.fqidx/.bai)
      --output-dir <DIR>          Directory to place all outputs (per-read.tsv, aggregate.tsv, etc.)
  -c, --contig <CONTIG>           Restrict processing to the provided contigs
      --contig-bins <TSV>         Optional TSV/CSV mapping contig identifiers to bins (contig<TAB>bin)
  -h, --help                      Print help
```