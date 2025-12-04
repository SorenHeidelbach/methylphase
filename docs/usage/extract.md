# Extract methylation calls

`methylphase extract` intersects motifs with mod-tagged alignments and emits per-read and aggregate methylation calls.

## Common arguments
- `--bam <mod.bam>`: indexed BAM with MM/ML tags.
- Motifs: `--motif CG_5,GATC_6mA_1` or `--motif-file motifs.tsv` (TSV/CSV with `motif`, `mod_type`, `mod_position`, optional `complement`, `bin`, `reference` columns).
- `--output-dir <dir>`: destination for TSVs (created if missing).
- `--sequence-fallback <fastq|bam>`: when SEQ is absent in BAM.
- Optional: `--motif-summary` to emit per-motif summaries; FASTQ subsets can be written when requested.

## Outputs
- `per_read.tsv`: methylation calls per read and motif occurrence.
- `aggregate.tsv`: aggregated motif-level calls.
- Optional `motif_summary.tsv` and FASTQ subsets if requested.

## Example
```bash
methylphase extract \
  --bam sample.mod.bam \
  --motif-file motifs.tsv \
  --output-dir results/extract \
  --motif-summary
```
