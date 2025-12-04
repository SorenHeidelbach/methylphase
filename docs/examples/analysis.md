# Example analysis

This walkthrough demonstrates a typical methylphase run on a mod-tagged BAM with motif-driven typing.

## Inputs
- `sample.mod.bam` + index
- Motifs file `motifs.tsv` (columns: `motif`, `mod_type`, `mod_position`, optional `complement`, `bin`, `reference`)
- Floria haploset `haplotypes.hapset`

## Steps
1) **Inspect contigs** (optional sanity check)
   ```bash
   methylphase utils contigs --bam sample.mod.bam
   ```
2) **Extract per-read methylation** (optional debugging artifact)
   ```bash
   methylphase extract \
     --bam sample.mod.bam \
     --motif-file motifs.tsv \
     --output-dir results/extract \
     --motif-summary
   ```
3) **Run the full phase-variants pipeline**
   ```bash
   methylphase phase-variants \
     --floria haplotypes.hapset \
     --bam sample.mod.bam \
     --motif-file motifs.tsv \
     --out results/phase_variants \
     --min-classes 1 \
     --max-classes 8
   ```
4) **Review clustering and model choice**
   - `results/phase_variants/split_reads/read_clustering_raw.tsv` for per-read motif features and clusters.
   - `results/phase_variants/fits/` for per-class log-likelihoods and criteria.
   - `results/phase_variants/summary.json` for the selected model and criteria values.
5) **Use the imputed outputs**
   - `imputed.tsv` for completed datasets (continuous and categorical fields).
   - `imputed_labels.tsv` for MAP class labels.
   - Feed these into downstream association tests or visualization tools.

## Optional follow-ups
- Convert aggregated methylation to VCF for external tools:
  ```bash
  methylphase utils vcf \
    --bam sample.mod.bam \
    --motif-file motifs.tsv \
    --output results/phase_variants/calls.vcf.bgz \
    --plain-output \
    --sample-name sample1
  tabix -p vcf results/phase_variants/calls.vcf.bgz
  ```
- Split FASTQs by inferred class for downstream assembly/analysis:
  ```bash
  methylphase typing split-fastq \
    --fastq reads.fastq.gz \
    --responsibilities results/phase_variants/best_responsibilities.tsv \
    --output-dir results/phase_variants/fastq_by_class
  ```
