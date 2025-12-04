# FAQ and troubleshooting

**The tool says MM/ML tags are missing.**  
Ensure your aligner wrote modification tags; some pipelines require `--mod-tags` or equivalent. Without MM/ML tags, methylphase cannot locate methylation coordinates.

**The BAM is missing an index.**  
Create one with `samtools index sample.mod.bam` (or provide a CRAM with a `.crai`).

**SEQ is absent in the BAM header.**  
Provide `--sequence-fallback <fastq|bam>` (and index the fallback if large) so methylphase can recover base sequences for motif placement.

**Clusters look noisy or over-split.**  
Try `--cluster-algorithm hdbscan` in `split-reads`, or tighten class bounds/penalty in `phase-variants` (e.g., `--max-classes 6`, `--penalty-multiplier 3.0`).

**Outputs are huge.**  
Disable FASTQ/BAM emission during `split-reads` unless needed (`--emit-fastq=false`), and limit motifs to those required for the analysis.

**Model selection feels off.**  
Inspect `fits/` and `summary.json` for criterion values; switch from `icl` to `bic` or `cv` (`--criterion`) and adjust `--penalty-multiplier` when models under/over fit.

**Performance is slow.**  
Increase `--threads` for extraction/clustering; keep BAM on fast storage and avoid overly large motif sets when possible.
