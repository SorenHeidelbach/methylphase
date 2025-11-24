import os
from glob import glob

configfile: "config.yaml"

SAMPLES = config.get("samples", {})
if not SAMPLES:
    raise ValueError("Add at least one sample under samples: in config.yaml")
DEFAULTS = config.get("defaults", {})
RESULTS_DIR = config.get("results_dir", "results")
ENV_CFG = config.get("env", {})
FLORIA_HOME = ENV_CFG.get("floria_bin", os.path.expanduser("~/.conda/envs/floria/bin"))
MODKIT_BIN = os.path.join(FLORIA_HOME, "modkit")
NANOMOTIF_BIN = os.path.join(FLORIA_HOME, "nanomotif")
FLORIA_BIN = os.path.join(FLORIA_HOME, "floria")


def cfg(sample, key, default_value):
    sample_cfg = SAMPLES.get(sample, {})
    if key in sample_cfg:
        return sample_cfg[key]
    if default_value is None:
        raise KeyError(f"No default for {key}; set samples.{sample}.{key} in config.yaml")
    value = default_value(sample) if callable(default_value) else default_value
    if isinstance(value, str):
        return value.format(sample=sample)
    return value


def result(sample, path):
    return os.path.join(RESULTS_DIR, sample, path)


def emit_fastq(sample):
    return SAMPLES.get(sample, {}).get(
        "split_emit_fastq", DEFAULTS.get("split_emit_fastq", False)
    )


def threads(sample, name, fallback):
    return SAMPLES.get(sample, {}).get(
        f"{name}_threads", DEFAULTS.get(f"{name}_threads", fallback)
    )


def _cluster_dir(sample):
    return os.path.join(result(sample, "split"), "clusters")


def _cluster_labels(sample):
    checkpoints.split_reads.get(sample=sample)
    pattern = os.path.join(_cluster_dir(sample), "cluster_*.bam")
    files = sorted(glob(pattern))
    if not files:
        raise ValueError(f"No cluster BAMs found for sample {sample}")
    return [os.path.splitext(os.path.basename(p))[0] for p in files]


def cluster_nanmotif_inputs(wildcards):
    clusters = _cluster_labels(wildcards.sample)
    return expand(
        "results/{sample}/clusters/{cluster}/nanomotif/bin-motifs.tsv",
        sample=wildcards.sample,
        cluster=clusters,
    )


FINAL_TARGETS = []
for sample in sorted(SAMPLES):
    FINAL_TARGETS.extend(
        [
            result(sample, "comparison/phasing_read_assignments.tsv"),
            result(sample, "comparison/phasing_pairwise_summary.tsv"),
            result(sample, "comparison/phasing_pairwise_optimal_mapping.tsv"),
            result(sample, "comparison/phasing_plots_similarity_heatmap.png"),
            result(sample, "comparison/phasing_plots_upset.png"),
            result(sample, "clusters/.nanomotif_complete"),
        ]
    )


rule all:
    input:
        "workflow.done"


rule workflow_done:
    input:
        FINAL_TARGETS
    output:
        touch("workflow.done")


rule modkit_pileup:
    input:
        bam=lambda w: cfg(w.sample, "bam", "data/{sample}/{sample}.bam"),
        reference=lambda w: cfg(w.sample, "reference", "data/{sample}/{sample}.fasta"),
    output:
        "results/{sample}/modkit/pileup.bed"
    params:
        threads=lambda w: threads(w.sample, "modkit", 4),
        extra=lambda w: SAMPLES.get(w.sample, {}).get("modkit_extra", DEFAULTS.get("modkit_extra", "")),
    shell:
        """
        mkdir -p $(dirname {output})
        {MODKIT_BIN} pileup "{input.bam}" "{output}" --reference "{input.reference}" \
            --threads {params.threads} {params.extra}
        """


rule nanomotif:
    input:
        assembly=lambda w: cfg(w.sample, "reference", "data/{sample}/{sample}.fasta"),
        pileup="results/{sample}/modkit/pileup.bed",
        contig_bins=lambda w: cfg(w.sample, "contig_bins", "config/contig_bins/{sample}.tsv"),
    output:
        "results/{sample}/nanomotif/bin-motifs.tsv"
    params:
        threads=lambda w: threads(w.sample, "nanomotif", 4),
        extra=lambda w: SAMPLES.get(w.sample, {}).get("nanomotif_extra", DEFAULTS.get("nanomotif_extra", ""))
    shell:
        """
        mkdir -p $(dirname {output})
        {NANOMOTIF_BIN} motif_discovery "{input.assembly}" "{input.pileup}" -c "{input.contig_bins}" \
            -t {params.threads} --out $(dirname {output}) {params.extra}
        """


rule extract:
    input:
        bam=lambda w: cfg(w.sample, "bam", "data/{sample}/{sample}.bam"),
        motif_file="results/{sample}/nanomotif/bin-motifs.tsv"
    output:
        per="results/{sample}/extract/per_read.tsv",
        agg="results/{sample}/extract/aggregate.tsv",
        flag="results/{sample}/extract/.extract_complete"
    params:
        fastq_dir=lambda w: result(w.sample, "extract/fastqs")
    shell:
        """
        mkdir -p {params.fastq_dir} $(dirname {output.per})
        methylation_phasing extract "{input.bam}" --motif-file "{input.motif_file}" \
            --per-read-tsv {output.per} --aggregate-tsv {output.agg} \
            --fastq-dir {params.fastq_dir}
        touch {output.flag}
        """


checkpoint split_reads:
    input:
        bam=lambda w: cfg(w.sample, "bam", "data/{sample}/{sample}.bam"),
        motif_file="results/{sample}/nanomotif/bin-motifs.tsv"
    output:
        clustering="results/{sample}/split/read_clustering.tsv",
        clustering_raw="results/{sample}/split/read_clustering_raw.tsv",
        summary="results/{sample}/split/cluster_summary.tsv",
        clusters_done="results/{sample}/split/.clusters_done"
    params:
        flag=lambda w: "--emit-fastq" if emit_fastq(w.sample) else "",
        outdir=lambda w: result(w.sample, "split")
    shell:
        """
        mkdir -p {params.outdir}
        methylation_phasing split-reads "{input.bam}" --output-dir {params.outdir} \
            --motif-file "{input.motif_file}" {params.flag}
        touch {output.clusters_done}
        """

rule cluster_get_fastq_with_mods:
    input:
        "results/{sample}/split/clusters/{cluster}.bam",
    output:
        "results/{sample}/clusters/{cluster}/{sample}_{cluster}.fastq"
    params:
        cluster=lambda w: w.cluster,
    shell:
        """
        mkdir -p $(dirname {output})
        samtools fastq -T MM,ML "{input}" > "{output}"
        """

rule cluster_assemble_myloasm:
    input:
        "results/{sample}/clusters/{cluster}/{sample}_{cluster}.fastq"
    output:
        directory=directory("results/{sample}/clusters/{cluster}/{sample}_{cluster}_myloasm"),
        asm="results/{sample}/clusters/{cluster}/{sample}_{cluster}_myloasm/assembly_primary.fa"
    shell:
        """
        myloasm "{input}" -o "{output.directory}" 
        """

rule cluster_contig_bin:
    input:
        assembly="results/{sample}/clusters/{cluster}/{sample}_{cluster}_myloasm/assembly_primary.fa"
    output:
        "results/{sample}/clusters/{cluster}/contig_bins.tsv"
    shell:
        """
        # take wildcard sample name as bin and assign all contigs to that bin
        mkdir -p $(dirname {output})
        samtools faidx "{input.assembly}"
        cut -f1 "{input.assembly}.fai" | awk -v bin="{wildcards.sample}_{wildcards.cluster}" '{{print $0"\t"bin}}' > "{output}"
        """

rule cluster_map:
    input:
        assembly="results/{sample}/clusters/{cluster}/{sample}_{cluster}_myloasm/assembly_primary.fa",
        fastq="results/{sample}/clusters/{cluster}/{sample}_{cluster}.fastq"
    output:
        bam="results/{sample}/clusters/{cluster}/{sample}_{cluster}.bam",
        bai="results/{sample}/clusters/{cluster}/{sample}_{cluster}.bam.bai"
    params:
        threads=lambda w: threads(w.sample, "modkit", 4),
    shell:
        """
        minimap2 -ayY -x map-ont -t {params.threads} "{input.assembly}" "{input.fastq}" | \
            samtools sort -o "{output.bam}"
        samtools index "{output.bam}"
        """

rule cluster_pileup:
    input:
        bam="results/{sample}/clusters/{cluster}/{sample}_{cluster}.bam",
        index="results/{sample}/clusters/{cluster}/{sample}_{cluster}.bam.bai"
    output:
        "results/{sample}/clusters/{cluster}/{sample}_{cluster}_pileup.bed"
    params:
        threads=lambda w: threads(w.sample, "modkit", 4),
        extra=lambda w: SAMPLES.get(w.sample, {}).get("modkit_extra", DEFAULTS.get("modkit_extra", "")),
    shell:
        """
        mkdir -p $(dirname {output})
        {MODKIT_BIN} pileup "{input.bam}" "{output}"  \
            --threads {params.threads} {params.extra}
        """

rule cluster_nanomotif:
    input:
        assembly="results/{sample}/clusters/{cluster}/{sample}_{cluster}_myloasm/assembly_primary.fa",
        pileup="results/{sample}/clusters/{cluster}/{sample}_{cluster}_pileup.bed",
        contig_bins="results/{sample}/clusters/{cluster}/contig_bins.tsv",
    output:
        "results/{sample}/clusters/{cluster}/nanomotif/bin-motifs.tsv"
    params:
        threads=lambda w: threads(w.sample, "nanomotif", 4),
        extra=lambda w: SAMPLES.get(w.sample, {}).get("nanomotif_extra", DEFAULTS.get("nanomotif_extra", "")),
    shell:
        """
        mkdir -p $(dirname {output})
        {NANOMOTIF_BIN} motif_discovery "{input.assembly}" "{input.pileup}" -c "{input.contig_bins}" \
            -t {params.threads} --out $(dirname {output}) {params.extra} -v
        """


rule cluster_nanmotif_complete:
    input:
        cluster_nanmotif_inputs
    output:
        "results/{sample}/clusters/.nanomotif_complete"
    shell:
        """
        mkdir -p $(dirname {output})
        touch {output}
        """


rule longitudinal:
    input:
        bam=lambda w: cfg(w.sample, "bam", "data/{sample}/{sample}.bam"),
        motif_file="results/{sample}/nanomotif/bin-motifs.tsv"
    output:
        "results/{sample}/longitudinal/longitudinal.tsv"
    shell:
        """
        mkdir -p $(dirname {output})
        methylation_phasing longitudinal "{input.bam}" --output-dir $(dirname {output}) \
            --motif-file "{input.motif_file}"
        """


rule vcf:
    input:
        bam=lambda w: cfg(w.sample, "bam", "data/{sample}/{sample}.bam"),
        motif_file="results/{sample}/nanomotif/bin-motifs.tsv"
    output:
        "results/{sample}/vcf/variants.vcf"
    shell:
        """
        mkdir -p $(dirname {output})
        methylation_phasing vcf "{input.bam}" --output {output} --motif-file "{input.motif_file}" --plain-output
        """


rule impute_bam:
    input:
        bam=lambda w: cfg(w.sample, "bam", "data/{sample}/{sample}.bam"),
        motif_file="results/{sample}/nanomotif/bin-motifs.tsv"
    output:
        bam="results/{sample}/impute/imputed.bam",
        summary="results/{sample}/impute/imputed.summary.tsv"
    params:
        threshold=lambda w: SAMPLES.get(w.sample, {}).get(
            "impute_threshold", DEFAULTS.get("impute_threshold", 0.8)
        )
    shell:
        """
        mkdir -p $(dirname {output.bam})
        methylation_phasing impute-bam "{input.bam}" --methylation-threshold {params.threshold} \
            --output {output.bam} --summary {output.summary} --motif-file "{input.motif_file}"
        samtools index {output.bam}
        """


rule longshot_imputed:
    input:
        bam="results/{sample}/impute/imputed.bam",
        reference=lambda w: cfg(w.sample, "reference", "data/{sample}/{sample}.fasta")
    output:
        "results/{sample}/longshot/imputed.vcf"
    shell:
        """
        mkdir -p $(dirname {output})
        samtools faidx "{input.reference}"
        longshot --bam "{input.bam}" --ref "{input.reference}" --out "{output}"
        """


rule longshot_original:
    input:
        bam=lambda w: cfg(w.sample, "bam", "data/{sample}/{sample}.bam"),
        reference=lambda w: cfg(w.sample, "reference", "data/{sample}/{sample}.fasta")
    output:
        "results/{sample}/longshot/original.vcf"
    shell:
        """
        mkdir -p $(dirname {output})
        samtools faidx "{input.reference}"
        longshot --bam "{input.bam}" --ref "{input.reference}" --out "{output}"
        """


rule floria_imputed:
    input:
        bam="results/{sample}/impute/imputed.bam",
        vcf="results/{sample}/longshot/imputed.vcf",
        reference=lambda w: cfg(w.sample, "reference", "data/{sample}/{sample}.fasta")
    output:
        "results/{sample}/floria_imputed/haploset.tsv"
    params:
        outdir=lambda w: result(w.sample, "floria_imputed/raw"),
        contig=lambda w: cfg(w.sample, "contig", "{sample}")
    shell:
        """
        mkdir -p {params.outdir}
        {FLORIA_BIN} -b "{input.bam}" -r "{input.reference}" -v "{input.vcf}" -o "{params.outdir}" --overwrite
        if [ -n "{params.contig}" ] && [ -f "{params.outdir}/{params.contig}/{params.contig}.haplosets" ]; then \
            cp "{params.outdir}/{params.contig}/{params.contig}.haplosets" "{output}"; \
        else \
            target=$(find "{params.outdir}" -name '*.haplosets' | head -n1) && cp "$target" "{output}"; \
        fi
        """


rule floria_original:
    input:
        bam=lambda w: cfg(w.sample, "bam", "data/{sample}/{sample}.bam"),
        vcf="results/{sample}/longshot/original.vcf",
        reference=lambda w: cfg(w.sample, "reference", "data/{sample}/{sample}.fasta")
    output:
        "results/{sample}/floria_original/haploset.tsv"
    params:
        outdir=lambda w: result(w.sample, "floria_original/raw"),
        contig=lambda w: cfg(w.sample, "contig", "{sample}")
    shell:
        """
        mkdir -p {params.outdir}
        {FLORIA_BIN} -b "{input.bam}" -r "{input.reference}" -v "{input.vcf}" -o "{params.outdir}" --overwrite
        if [ -n "{params.contig}" ] && [ -f "{params.outdir}/{params.contig}/{params.contig}.haplosets" ]; then \
            cp "{params.outdir}/{params.contig}/{params.contig}.haplosets" "{output}"; \
        else \
            target=$(find "{params.outdir}" -name '*.haplosets' | head -n1) and cp "$target" "{output}"; \
        fi
        """



rule compare_phasing:
    input:
        floria_imputed="results/{sample}/floria_imputed/haploset.tsv",
        floria_original="results/{sample}/floria_original/haploset.tsv",
        split="results/{sample}/split/read_clustering.tsv"
    output:
        assignments="results/{sample}/comparison/phasing_read_assignments.tsv",
        summary="results/{sample}/comparison/phasing_pairwise_summary.tsv",
        optimal="results/{sample}/comparison/phasing_pairwise_optimal_mapping.tsv",
        heatmap="results/{sample}/comparison/phasing_plots_similarity_heatmap.png",
        upset="results/{sample}/comparison/phasing_plots_upset.png"
    params:
        prefix=lambda w: result(w.sample, "comparison/phasing")
    shell:
        """
        mkdir -p $(dirname {output.assignments})
        Rscript scripts/compare_phasing.R --output-prefix "{params.prefix}" \
            --input imputed-floria={input.floria_imputed} \
            --input floria={input.floria_original} \
            --input meth-cluster={input.split}
        """
