mod binning;
pub mod contigs;
pub mod extract;
pub mod impute_bam;
pub mod longitudinal;
pub mod reads;
pub mod shared;
pub mod split_reads;
pub mod vcf;

use crate::cli::{Cli, Command, ContigSelection};
use anyhow::Result;
use std::{
    collections::HashSet,
    path::{Path, PathBuf},
};

use self::binning::{BinAssignment, ContigBins};

pub fn run(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Contigs { bam, contig_args } => {
            run_with_bins("contigs", contig_args, move |_bin, contigs| {
                let bam = bam.clone();
                contigs::run(bam, contigs)
            })
        }
        Command::Reads {
            bam,
            contig_args,
            limit,
        } => run_with_bins("reads", contig_args, move |_bin, contigs| {
            let bam = bam.clone();
            reads::run(bam, contigs, limit)
        }),
        Command::Extract {
            bam,
            motifs,
            motif_file,
            motif_summary_tsv,
            fastq_dir,
            sequence_fallback,
            sequence_index,
            output_dir,
            contig_args,
        } => run_with_bins("extract", contig_args, move |bin, contigs| {
            let bam = bam.clone();
            let motifs = motifs.clone();
            let motif_file = motif_file.clone();
            let motif_summary_tsv = motif_summary_tsv.clone();
            let fastq_dir = fastq_dir.clone();
            let sequence_fallback = sequence_fallback.clone();
            let sequence_index = sequence_index.clone();
            let mut output_dir = output_dir.clone();
            if let Some(bin) = bin {
                output_dir = apply_bin_directory(&output_dir, bin);
            }
            let bin_ids = bin.map(bin_identifier_values);
            extract::run(
                bam,
                motifs,
                motif_file,
                motif_summary_tsv,
                fastq_dir,
                sequence_fallback,
                sequence_index,
                output_dir,
                bin_ids,
                contigs,
            )
        }),
        Command::Longitudinal {
            bam,
            motifs,
            motif_file,
            block_size,
            methylation_threshold,
            sequence_fallback,
            sequence_index,
            output_dir,
            contig_args,
        } => run_with_bins("longitudinal", contig_args, move |bin, contigs| {
            let bam = bam.clone();
            let motifs = motifs.clone();
            let motif_file = motif_file.clone();
            let sequence_fallback = sequence_fallback.clone();
            let sequence_index = sequence_index.clone();
            let mut output_dir = output_dir.clone();
            if let Some(bin) = bin {
                output_dir = apply_bin_directory(&output_dir, bin);
            }
            let bin_ids = bin.map(bin_identifier_values);
            longitudinal::run(
                bam,
                motifs,
                motif_file,
                block_size,
                methylation_threshold,
                sequence_fallback,
                sequence_index,
                output_dir,
                bin_ids,
                contigs,
            )
        }),
        Command::SplitReads {
            bam,
            motifs,
            motif_file,
            sequence_fallback,
            sequence_index,
            output_dir,
            min_cluster_size,
            min_samples,
            emit_fastq,
            threads,
            cluster_algorithm,
            contig_args,
        } => run_with_bins("split_reads", contig_args, move |bin, contigs| {
            let bam = bam.clone();
            let motifs = motifs.clone();
            let motif_file = motif_file.clone();
            let sequence_fallback = sequence_fallback.clone();
            let sequence_index = sequence_index.clone();
            let mut output_dir = output_dir.clone();
            if let Some(bin) = bin {
                output_dir = apply_bin_directory(&output_dir, bin);
            }
            let bin_ids = bin.map(bin_identifier_values);
            split_reads::run(
                bam,
                motifs,
                motif_file,
                sequence_fallback,
                sequence_index,
                output_dir,
                min_cluster_size,
                min_samples,
                emit_fastq,
                threads,
                cluster_algorithm,
                bin_ids,
                contigs,
            )
        }),
        Command::Vcf {
            bam,
            motifs,
            motif_file,
            sequence_fallback,
            sequence_index,
            output,
            sample_name,
            methylation_threshold,
            plain_output,
            contig_args,
        } => run_with_bins("vcf", contig_args, move |bin, contigs| {
            let bam = bam.clone();
            let motifs = motifs.clone();
            let motif_file = motif_file.clone();
            let sequence_fallback = sequence_fallback.clone();
            let sequence_index = sequence_index.clone();
            let mut output = output.clone();
            if let Some(bin) = bin {
                output = apply_bin_file(&output, bin);
            }
            let bin_ids = bin.map(bin_identifier_values);
            vcf::run(
                bam,
                motifs,
                motif_file,
                sequence_fallback,
                sequence_index,
                output,
                sample_name.clone(),
                methylation_threshold,
                plain_output,
                bin_ids,
                contigs,
            )
        }),
        Command::ImputeBam {
            bam,
            output,
            methylation_threshold,
            summary,
            motifs,
            motif_file,
            impute_all,
            contig_args,
        } => run_with_bins("impute_bam", contig_args, move |bin, contigs| {
            let bam = bam.clone();
            let motifs = motifs.clone();
            let motif_file = motif_file.clone();
            let mut output_path = output.clone();
            if let Some(bin) = bin {
                output_path = apply_bin_file(&output_path, bin);
            }
            let summary_path = summary.as_ref().map(|path| match bin {
                Some(bin) => apply_bin_file(path, bin),
                None => path.clone(),
            });
            let bin_ids = bin.map(bin_identifier_values);
            impute_bam::run(
                bam,
                output_path,
                methylation_threshold,
                summary_path,
                motifs,
                motif_file,
                impute_all,
                bin_ids,
                contigs,
            )
        }),
    }
}

fn bin_identifier_values(bin: &BinAssignment) -> Vec<String> {
    let mut values = vec![bin.name().to_string()];
    let safe = bin.safe_name();
    if !safe.eq(bin.name()) {
        values.push(safe.to_string());
    }
    values
}

fn run_with_bins<F>(command: &str, selection: ContigSelection, mut op: F) -> Result<()>
where
    F: FnMut(Option<&BinAssignment>, Vec<String>) -> Result<()>,
{
    let ContigSelection {
        contigs,
        contig_bins,
    } = selection;

    if let Some(path) = contig_bins {
        let bins = ContigBins::from_path(&path)?;
        for bin in bins.iter() {
            if let Some(selected) = prepare_contigs(command, Some(bin), &contigs) {
                op(Some(bin), selected)?;
            }
        }
        Ok(())
    } else if let Some(selected) = prepare_contigs(command, None, &contigs) {
        op(None, selected)
    } else {
        Ok(())
    }
}

fn prepare_contigs(
    command: &str,
    bin: Option<&BinAssignment>,
    requested: &[String],
) -> Option<Vec<String>> {
    match contigs_for_execution(bin, requested) {
        Some(values) => {
            if let Some(bin) = bin {
                eprintln!(
                    "{command}: processing bin {} with {} contigs",
                    bin.name(),
                    values.len()
                );
            }
            Some(values)
        }
        None => {
            if let Some(bin) = bin {
                eprintln!(
                    "{command}: skipping bin {} (no matching contigs)",
                    bin.name()
                );
            }
            None
        }
    }
}

fn contigs_for_execution(bin: Option<&BinAssignment>, requested: &[String]) -> Option<Vec<String>> {
    match bin {
        Some(bin) => {
            if requested.is_empty() {
                if bin.contigs().is_empty() {
                    None
                } else {
                    Some(bin.contigs().to_vec())
                }
            } else {
                let allowed: HashSet<&str> = bin.contigs().iter().map(|c| c.as_str()).collect();
                let filtered: Vec<String> = requested
                    .iter()
                    .filter(|c| allowed.contains(c.as_str()))
                    .cloned()
                    .collect();
                if filtered.is_empty() {
                    None
                } else {
                    Some(filtered)
                }
            }
        }
        None => Some(requested.to_vec()),
    }
}

fn apply_bin_directory(base: &PathBuf, bin: &BinAssignment) -> PathBuf {
    base.join(bin.safe_name())
}

fn apply_bin_file(base: &PathBuf, bin: &BinAssignment) -> PathBuf {
    append_suffix(base.as_path(), bin.safe_name())
}

fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let parent = path
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(PathBuf::new);
    let file_name = path
        .file_name()
        .map(|name| name.to_string_lossy().to_string())
        .unwrap_or_default();

    let new_name = if file_name.is_empty() {
        suffix.to_string()
    } else if let Some((stem, ext)) = file_name.rsplit_once('.') {
        if stem.is_empty() {
            format!("{}.{}", suffix, ext)
        } else {
            format!("{stem}.{suffix}.{ext}")
        }
    } else {
        format!("{file_name}.{suffix}")
    };

    parent.join(new_name)
}
