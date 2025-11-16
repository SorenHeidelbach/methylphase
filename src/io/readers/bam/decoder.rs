use crate::core::{ReadRecord, Strand};
use crate::io::readers::bam::modifications;
use anyhow::{Context, Result};
use bstr::ByteSlice;
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles_sam::{self as sam};
use std::cell::Cell;

const MAX_WARNINGS: usize = 20;

/// Converts BAM records into `ReadRecord` instances without requiring the reader
/// to be tied to the downstream data model.
pub struct BamRecordDecoder {
    reference_names: Vec<String>,
    parse_failures: Cell<usize>,
}

impl BamRecordDecoder {
    pub fn new(header: &sam::Header) -> Self {
        let reference_names = header
            .reference_sequences()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect();

        Self {
            reference_names,
            parse_failures: Cell::new(0),
        }
    }

    pub fn decode(&self, record: bam::Record) -> Result<ReadRecord> {
        let name_field = record.name().context("encountered read without name")?;
        let read_id = name_field
            .to_str()
            .context("read name contains invalid UTF-8")?
            .to_string();

        let contig = record
            .reference_sequence_id()
            .transpose()?
            .and_then(|idx| self.reference_names.get(idx))
            .cloned();

        let flags = record.flags();
        let strand = if flags.is_unmapped() {
            None
        } else if flags.is_reverse_complemented() {
            Some(Strand::Reverse)
        } else {
            Some(Strand::Forward)
        };

        let sequence: Vec<u8> = record.sequence().iter().collect();
        let quality_scores = record.quality_scores().as_ref().to_vec();
        let reference_positions = map_read_to_reference(&record)?;
        let modifications = match modifications::extract(&record) {
            Ok(calls) => calls,
            Err(err) => {
                self.record_parse_failure(&read_id, err);
                Vec::new()
            }
        };

        Ok(ReadRecord::new(
            read_id,
            contig,
            strand,
            sequence,
            modifications,
            quality_scores,
            reference_positions,
        ))
    }

    pub fn report(&self) {
        let skipped = self.parse_failures.get();
        if skipped > 0 {
            eprintln!(
                "warning: skipped {} reads with incompatible MM/ML tags (counts suppressed after first {}).",
                skipped,
                MAX_WARNINGS
            );
        }
    }

    fn record_parse_failure(&self, read_id: &str, err: anyhow::Error) {
        let count = self.parse_failures.get();
        if count < MAX_WARNINGS {
            eprintln!(
                "warning: failed to parse modifications for read {}: {}",
                read_id, err
            );
        } else if count == MAX_WARNINGS {
            eprintln!("warning: additional modification parse failures suppressed");
        }
        self.parse_failures.set(count.saturating_add(1));
    }
}

fn map_read_to_reference(record: &bam::Record) -> Result<Vec<Option<i64>>> {
    let read_len = record.sequence().len();
    let mut mapping = vec![None; read_len];

    let Some(start) = record.alignment_start().transpose()? else {
        return Ok(mapping);
    };

    let mut read_idx = 0usize;
    let mut ref_pos = start.get() as i64;

    for op_result in record.cigar().iter() {
        let op = op_result?;
        let len = op.len() as usize;
        match op.kind() {
            CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                for _ in 0..len {
                    if read_idx < mapping.len() {
                        mapping[read_idx] = Some(ref_pos);
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
