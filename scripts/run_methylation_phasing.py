#!/usr/bin/env python3
"""
Wrapper CLI to impute BAM methylation, emit a VCF, and run Whatshap phasing.
"""

import argparse
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List


def run_cmd(cmd: List[str]) -> None:
    """Execute a command, printing it beforehand."""
    print("+", " ".join(shlex.quote(part) for part in cmd))
    subprocess.run(cmd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run methylphase impute/VCF steps followed by Whatshap polyphase."
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="Input BAM with MM/ML tags (will remain unchanged).",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA passed to Whatshap.",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for generated files (imputed BAM, VCF, phased VCF).",
    )
    parser.add_argument(
        "--motif",
        dest="motifs",
        action="append",
        default=[],
        help=(
            "Motif descriptor(s) (motif_modtype_modpos). Repeat flag or use commas "
            "to supply multiple motifs."
        ),
    )
    parser.add_argument(
        "--motif-file",
        help="Optional TSV/CSV describing motifs (same format as methylphase).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.8,
        help="Methylation probability threshold for imputing bases (default: 0.8).",
    )
    parser.add_argument(
        "--sample-name",
        help="Sample name for the VCF header (defaults to output prefix basename).",
    )
    parser.add_argument(
        "--plain-vcf",
        action="store_true",
        help="Write an uncompressed VCF (bgzip/tabix skipped).",
    )
    parser.add_argument(
        "--whatshap-bin",
        default="whatshap",
        help="Path to the Whatshap executable (default: whatshap on PATH).",
    )
    parser.add_argument(
        "--whatshap-extra",
        nargs=argparse.REMAINDER,
        default=[],
        help="Additional arguments passed to Whatshap (use after `--`).",
    )
    parser.add_argument(
        "--ploidy",
        type=int,
        required=True,
        help="Ploidy to pass to whatshap polyphase (--ploidy).",
    )

    args = parser.parse_args()

    if not args.motifs and not args.motif_file:
        parser.error("must supply at least one --motif (or comma list) or a --motif-file")

    args.combined_motifs = ",".join(args.motifs)
    return args


def main() -> int:
    args = parse_args()

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    imputed_bam = prefix.with_suffix(".imputed.bam")
    vcf_suffix = ".vcf" if args.plain_vcf else ".vcf.gz"
    vcf_path = prefix.with_suffix(vcf_suffix)
    phased_vcf = prefix.with_suffix(".phased.vcf")

    sample_name = args.sample_name or prefix.name

    # Step 1: impute high-confidence methylation into the BAM.
    impute_cmd = [
        "methylphase",
        "impute-bam",
        "--methylation-threshold",
        str(args.threshold),
        "--output",
        str(imputed_bam),
        args.bam,
    ]
    run_cmd(impute_cmd)

    # Step 2: emit methylation VCF.
    vcf_cmd = [
        "methylphase",
        "vcf",
        "--output",
        str(vcf_path),
        "--sample-name",
        sample_name,
        "--methylation-threshold",
        str(args.threshold),
    ]
    if args.plain_vcf:
        vcf_cmd.append("--plain-output")
    if args.combined_motifs:
        vcf_cmd.extend(["--motif", args.combined_motifs])
    if args.motif_file:
        vcf_cmd.extend(["--motif-file", args.motif_file])
    vcf_cmd.append("--")
    vcf_cmd.append(args.bam)
    run_cmd(vcf_cmd)

    # Step 3: run Whatshap polyphase with methylation VCF + imputed BAM.
    whatshap_cmd = [
        args.whatshap_bin,
        "polyphase",
        "--reference",
        args.reference,
        "--ploidy",
        str(args.ploidy),
        "-o",
        str(phased_vcf),
    ]
    if args.whatshap_extra:
        whatshap_cmd.extend(args.whatshap_extra)
    whatshap_cmd.extend([str(vcf_path), str(imputed_bam)])
    run_cmd(whatshap_cmd)

    print("Pipeline complete.")
    print(f"- Imputed BAM: {imputed_bam}")
    print(f"- Methylation VCF: {vcf_path}")
    print(f"- Whatshap output: {phased_vcf}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as err:
        cmd = " ".join(shlex.quote(part) for part in err.cmd)
        print(f"Command failed ({err.returncode}): {cmd}", file=sys.stderr)
        raise
