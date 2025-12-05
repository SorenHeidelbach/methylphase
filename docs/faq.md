---
layout: default
title: FAQ and troubleshooting
---

# FAQ and troubleshooting

**How do i generate Floria haplotypes?**  
Use the [Floria toolkit](https://github.com/bluenote-1577/floria) to call haplotypes from variant VCFs and BAMs. See the [Floria docs](https://phase-doc.readthedocs.io/en/latest/introduction.html) for details. 

**The tool says MM/ML tags are missing.**  
Ensure your aligner wrote modification tags; some pipelines require `--mod-tags` or equivalent. Without MM/ML tags, methylphase cannot locate methylation coordinates. We recommend mapping with [dorado aligner](https://software-docs.nanoporetech.com/dorado/latest/basecaller/alignment/) using BAM files with MM/ML tags. If you prefer minimap2, use `-y` flag to ensure read-id flags are emitted to the resulting alignment BAM.

**The BAM is missing an index.**  
Create one with `samtools index sample.mod.bam`.

**SEQ is absent in the BAM header.**  
Provide `--sequence-fallback <fastq|bam>` (and index the fallback if large) so methylphase can recover base sequences for motif placement.

**Too many or few variant calls.**  
If you observe oversplitting into variants, increase penalty in `phase-variants` (e.g., `--penalty-multiplier 3.0`). If you observe too few variants, consider changing the default criterion. The ICL criterion tends produce the least amount of variants; switching to BIC or CV (`--criterion bic` or `--criterion cv`) may increase variant calls.

**Outputs are huge.**  
Disable FASTQ/BAM emission during `split-reads` unless needed (`--emit-fastq=false`), and limit motifs to those required for the analysis.

**How to investigate variant selection.**  
Inspect `fits/` and `summary.json` for Pi values, these indicate the relative abundance of the variants. Low Pi variants (e.g., < 0.05) may be spurious and can be filtered post-hoc or rerun with adjusted parameters.
