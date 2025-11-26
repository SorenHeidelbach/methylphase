use crate::{
    commands::shared,
    core::{ModificationCall, MotifQuery},
    io::readers::bam::modifications,
};
use anyhow::{bail, Context, Result};
use bstr::ByteSlice;
use noodles_bam as bam;
use noodles_sam::alignment::{
    io::Write as SamWrite,
    record::{cigar::op::Kind as CigarKind, Cigar as _},
    RecordBuf,
};
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

pub fn run(
    bam: PathBuf,
    output: PathBuf,
    methylation_threshold: f32,
    summary: Option<PathBuf>,
    motifs: Vec<String>,
    motif_file: Option<PathBuf>,
    impute_all: bool,
    bin_ids: Option<Vec<String>>,
    contigs: Vec<String>,
) -> Result<()> {
    if !(0.0..=1.0).contains(&methylation_threshold) {
        bail!("methylation-threshold must be between 0.0 and 1.0");
    }

    let bin_id_refs = bin_ids.as_deref();
    let motif_queries = if impute_all {
        None
    } else {
        let queries = shared::load_motif_queries(motifs, motif_file, bin_id_refs)?
            .into_iter()
            .collect::<Vec<_>>();
        if queries.is_empty() {
            bail!("no motifs provided; use --motif/--motif-file or --impute-all");
        }
        Some(queries)
    };

    let summary_path = summary.unwrap_or_else(|| default_summary_path(&output));
    let column_counts = build_column_counts(&bam)?;

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&bam)
        .with_context(|| format!("failed to open BAM {}", bam.display()))?;
    let header = reader.read_header().context("failed to read BAM header")?;

    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    writer
        .write_header(&header)
        .context("failed to write BAM header")?;

    let mut summary_writer = BufWriter::new(
        File::create(&summary_path)
            .with_context(|| format!("failed to create {}", summary_path.display()))?,
    );
    writeln!(
        summary_writer,
        "contig\treference_pos\tread_id\told_base\tnew_base\tprobability"
    )?;

    let reference_names: Vec<_> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    let contig_filter = if contigs.is_empty() {
        None
    } else {
        Some(contigs.into_iter().collect::<HashSet<_>>())
    };

    let mut record = RecordBuf::default();
    let mut processed_reads = 0usize;
    let mut imputed_calls = 0usize;
    let mut skipped_reads = 0usize;

    loop {
        let bytes = reader
            .read_record_buf(&header, &mut record)
            .context("failed to read BAM record")?;
        if bytes == 0 {
            break;
        }

        processed_reads += 1;
        let read_name = record
            .name()
            .map(|n| n.to_str_lossy().into_owned())
            .unwrap_or_else(|| "<unnamed>".to_string());

        let calls = match modifications::extract_from_record_buf(&record) {
            Ok(calls) => calls,
            Err(err) => {
                eprintln!(
                    "warning: skipping read {} due to MM/ML parse error: {}",
                    read_name, err
                );
                skipped_reads += 1;
                continue;
            }
        };

        let Some(ref_id) = record.reference_sequence_id() else {
            continue;
        };

        let contig_name = reference_names
            .get(ref_id)
            .cloned()
            .unwrap_or_else(|| "<unknown>".to_string());

        if let Some(filter) = contig_filter.as_ref() {
            if !filter.contains(&contig_name) {
                continue;
            }
        }

        let reference_positions = map_read_to_reference_positions(&record)?;
        let motif_allowed_positions = motif_queries
            .as_ref()
            .map(|motifs| build_allowed_positions(record.sequence().as_ref(), motifs));
        let mut replacements = 0usize;

        let quality_len = record.quality_scores().as_ref().len();
        let mut pending_quality = Vec::new();

        {
            let sequence = record.sequence_mut().as_mut();
            for call in calls {
                if let Some(allowed) = motif_allowed_positions.as_ref() {
                    match allowed.get(&call.position) {
                        Some(labels) if labels.contains(&call.kind.label) => {}
                        _ => continue,
                    }
                }

                if !should_impute(&call, methylation_threshold) {
                    continue;
                }

                if call.position >= sequence.len() || call.position >= reference_positions.len() {
                    continue;
                }

                let Some(reference_position) = reference_positions[call.position] else {
                    continue;
                };

                let key = ColumnKey {
                    ref_id,
                    position: reference_position,
                };

                let Some(column) = column_counts.get(&key) else {
                    continue;
                };

                let original_base = sequence[call.position].to_ascii_uppercase();
                if !is_standard_base(original_base) {
                    continue;
                }

                let Some(new_base) = choose_replacement(&column.counts, original_base) else {
                    continue;
                };

                if new_base == original_base {
                    continue;
                }

                sequence[call.position] = new_base;
                replacements += 1;

                if let Some(probability) = call.probability {
                    if let Some(score) = probability_to_quality(probability) {
                        if call.position < quality_len {
                            pending_quality.push((call.position, score));
                        }
                    }
                }

                writeln!(
                    summary_writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    contig_name,
                    reference_position,
                    read_name,
                    original_base as char,
                    new_base as char,
                    call.probability
                        .map(|p| format!("{:.4}", p))
                        .unwrap_or_else(|| "NA".to_string())
                )?;
            }
        }

        if replacements > 0 {
            imputed_calls += replacements;
            let qualities = record.quality_scores_mut().as_mut();
            for (position, score) in pending_quality {
                if position < qualities.len() {
                    qualities[position] = score;
                }
            }
        }

        writer
            .write_alignment_record(&header, &record)
            .context("failed to write BAM record")?;
    }

    summary_writer.flush()?;
    writer.try_finish().context("failed to finish BAM writer")?;

    eprintln!(
        "impute-bam: processed {} reads ({} skipped); imputed {} modification sites into {}",
        processed_reads,
        skipped_reads,
        imputed_calls,
        output.display()
    );
    eprintln!("impute-bam: summary written to {}", summary_path.display());

    Ok(())
}

fn default_summary_path(output: &Path) -> PathBuf {
    output.with_extension("summary.tsv")
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
struct ColumnKey {
    ref_id: usize,
    position: i64,
}

struct ColumnCounts {
    counts: [usize; 4],
}

fn build_column_counts(path: &Path) -> Result<HashMap<ColumnKey, ColumnCounts>> {
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("failed to open BAM {}", path.display()))?;
    let header = reader.read_header().context("failed to read BAM header")?;
    let mut record = RecordBuf::default();
    let mut counts: HashMap<ColumnKey, ColumnCounts> = HashMap::new();

    loop {
        let bytes = reader
            .read_record_buf(&header, &mut record)
            .context("failed to read BAM record")?;
        if bytes == 0 {
            break;
        }

        let Some(ref_id) = record.reference_sequence_id() else {
            continue;
        };

        let positions = map_read_to_reference_positions(&record)?;
        let sequence = record.sequence().as_ref();

        for (idx, base) in sequence.iter().enumerate() {
            let Some(position) = positions.get(idx).and_then(|p| *p) else {
                continue;
            };

            if let Some(base_idx) = base_index(*base) {
                let entry = counts
                    .entry(ColumnKey { ref_id, position })
                    .or_insert_with(|| ColumnCounts { counts: [0; 4] });
                entry.counts[base_idx] += 1;
            }
        }
    }

    Ok(counts)
}

fn build_allowed_positions(
    sequence: &[u8],
    motifs: &[MotifQuery],
) -> HashMap<usize, HashSet<String>> {
    let mut map: HashMap<usize, HashSet<String>> = HashMap::new();
    for motif in motifs {
        let offset = motif.mod_position();
        for start in motif.matches(sequence) {
            let pos = start + offset;
            map.entry(pos)
                .or_insert_with(HashSet::new)
                .insert(motif.mod_label().to_string());
        }
    }
    map
}

fn map_read_to_reference_positions(record: &RecordBuf) -> Result<Vec<Option<i64>>> {
    let mut mapping = vec![None; record.sequence().len()];
    let Some(start) = record.alignment_start() else {
        return Ok(mapping);
    };
    let mut read_idx = 0usize;
    let mut ref_pos = start.get() as i64;

    for op in record.cigar().iter() {
        let op = op?;
        let len = op.len() as usize;
        match op.kind() {
            CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                for _ in 0..len {
                    if let Some(slot) = mapping.get_mut(read_idx) {
                        *slot = Some(ref_pos);
                    }
                    read_idx += 1;
                    ref_pos += 1;
                }
            }
            CigarKind::Insertion | CigarKind::SoftClip => {
                read_idx += len;
            }
            CigarKind::Deletion | CigarKind::Skip => {
                ref_pos += len as i64;
            }
            CigarKind::HardClip | CigarKind::Pad => {}
        }
    }

    Ok(mapping)
}

fn should_impute(call: &ModificationCall, threshold: f32) -> bool {
    call.methylated && call.probability.map(|p| p >= threshold).unwrap_or(true)
}

fn choose_replacement(counts: &[usize; 4], avoid: u8) -> Option<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    for (idx, base) in bases.iter().enumerate() {
        if *base == avoid {
            continue;
        }
        if counts[idx] == 0 {
            return Some(*base);
        }
    }

    let mut best_idx = None;
    let mut best_count = usize::MAX;
    for (idx, count) in counts.iter().enumerate() {
        let base = bases[idx];
        if base == avoid {
            continue;
        }

        if *count < best_count {
            best_idx = Some(idx);
            best_count = *count;
        }
    }

    best_idx.map(|idx| bases[idx])
}

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn is_standard_base(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn probability_to_quality(prob: f32) -> Option<u8> {
    if !(0.0..=1.0).contains(&prob) {
        return None;
    }

    if (1.0 - prob).abs() < f32::EPSILON {
        return Some(93);
    }

    let clamped = prob.clamp(0.0, 1.0);
    let complement = (1.0 - clamped).max(1e-9);
    let phred = -10.0 * complement.log10();
    let rounded = phred.round().clamp(0.0, 93.0);
    Some(rounded as u8)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn probability_to_score_clamps() {
        let score = probability_to_quality(0.9).unwrap();
        assert!(score > 0);
        let high = probability_to_quality(1.0).unwrap();
        assert_eq!(high, 93);
        let zero = probability_to_quality(0.0).unwrap();
        assert_eq!(zero, 0);
    }

    #[test]
    fn choose_replacement_prefers_missing_base() {
        let counts = [5, 0, 0, 1];
        assert_eq!(choose_replacement(&counts, b'A'), Some(b'C'));
    }

    #[test]
    fn choose_replacement_uses_least_frequent_when_all_present() {
        let counts = [5, 3, 4, 6];
        assert_eq!(choose_replacement(&counts, b'A'), Some(b'C'));
    }
}
