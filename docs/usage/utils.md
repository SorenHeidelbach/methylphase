---
layout: default
title: Utility commands
---

# Utility commands

`methylphase utils <subcommand>` provides supporting tools for inspection and downstream use.

## vcf
Summarize motif methylation to VCF (optionally bgzip+tabix).
```bash
methylphase utils vcf \
  --bam sample.mod.bam \
  --motif CG_5 \
  --output calls.vcf.bgz \
  --plain-output \
  --sample-name sample1
```

## impute-bam
Impute bases/qualities directly into a BAM using methylation probabilities.
```bash
methylphase utils impute-bam \
  --bam sample.mod.bam \
  --motif CG_5 \
  --output imputed.bam \
  --prob-threshold 0.5
```

## contigs
List contigs in a BAM (with optional lengths).
```bash
methylphase utils contigs --bam sample.mod.bam
```

## reads
List read IDs for selected contigs (optionally capped).
```bash
methylphase utils reads --bam sample.mod.bam --contig chr1 --limit 1000
```

## longitudinal
Compute longitudinal methylation agreement across contigs/bins.
```bash
methylphase utils longitudinal \
  --bam sample.mod.bam \
  --motif-file motifs.tsv \
  --output longitudinal.tsv
```
