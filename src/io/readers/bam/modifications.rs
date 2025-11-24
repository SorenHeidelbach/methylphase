use crate::core::{ModificationCall, ModificationKind, Strand};
use anyhow::{anyhow, bail, Context, Result};
use bstr::ByteSlice;
use noodles_bam as bam;
use noodles_sam::alignment::{
    record::data::field::{value::Array, Tag, Value},
    record_buf::{self, data::field::value::Array as BufArray},
    RecordBuf,
};

const MM_TAG: Tag = Tag::BASE_MODIFICATIONS;
const ML_TAG: Tag = Tag::BASE_MODIFICATION_PROBABILITIES;

pub fn extract(record: &bam::Record) -> Result<Vec<ModificationCall>> {
    let data = record.data();
    let Some(mm_value) = data.get(&MM_TAG).transpose()? else {
        return Ok(Vec::new());
    };

    let mm = match mm_value {
        Value::String(s) => s.to_str().context("MM tag contains invalid UTF-8")?,
        _ => bail!("MM tag must be a string"),
    };

    let ml = parse_ml_field(data.get(&ML_TAG))?;
    let sequence: Vec<u8> = record.sequence().iter().collect();
    parse_modifications(mm, &sequence, ml.as_deref())
}

pub fn extract_from_record_buf(record: &RecordBuf) -> Result<Vec<ModificationCall>> {
    use record_buf::data::field::Value as BufValue;

    let data = record.data();
    let Some(mm_value) = data.get(&MM_TAG) else {
        return Ok(Vec::new());
    };

    let mm_owned = match mm_value {
        BufValue::String(s) => s
            .to_str()
            .map_err(|e| anyhow!("MM tag contains invalid UTF-8: {e}"))?
            .to_string(),
        _ => bail!("MM tag must be a string"),
    };

    let ml = match data.get(&ML_TAG) {
        Some(BufValue::Array(BufArray::UInt8(values))) => Some(values.clone()),
        Some(other) => bail!(
            "ML tag must be an array of uint8 values, got {:?}",
            other.ty()
        ),
        None => None,
    };

    let sequence = record.sequence().as_ref().to_vec();
    parse_modifications(&mm_owned, &sequence, ml.as_deref())
}

fn parse_ml_field(field: Option<Result<Value<'_>, std::io::Error>>) -> Result<Option<Vec<u8>>> {
    let Some(value) = field.transpose()? else {
        return Ok(None);
    };

    match value {
        Value::Array(Array::UInt8(values)) => {
            let probs = values
                .iter()
                .collect::<std::io::Result<Vec<_>>>()
                .context("invalid ML tag")?;
            Ok(Some(probs))
        }
        other => bail!(
            "ML tag must be an array of uint8 values, got {:?}",
            other.ty()
        ),
    }
}

fn parse_modifications(
    mm: &str,
    sequence: &[u8],
    ml: Option<&[u8]>,
) -> Result<Vec<ModificationCall>> {
    let segments = parse_segments(mm)?;
    let mut calls = Vec::new();
    let mut probability_index = 0;

    for segment in segments {
        let positions = canonical_positions(sequence, segment.canonical_base, segment.reverse);
        let mut cursor: isize = -1;
        let mut used_positions = vec![false; positions.len()];

        for delta in segment.deltas {
            cursor += delta as isize + 1;
            if cursor < 0 {
                bail!("negative modification offset encountered");
            }

            let canonical_idx = cursor as usize;
            let Some(read_pos) = positions.get(canonical_idx).copied() else {
                let base = char::from(segment.canonical_base);
                bail!(
                    "MM segment for base {} references canonical index {} but only {} occurrences were found",
                    base,
                    canonical_idx,
                    positions.len()
                );
            };

            let probability = if let Some(values) = ml {
                let value = values.get(probability_index).copied().ok_or_else(|| {
                    anyhow!(
                        "ML tag shorter than MM events (needed {})",
                        probability_index + 1
                    )
                })?;
                probability_index += 1;
                Some(value as f32 / 255.0)
            } else {
                None
            };

            if let Some(flag) = used_positions.get_mut(canonical_idx) {
                *flag = true;
            }

            calls.push(ModificationCall {
                position: read_pos,
                probability,
                methylated: true,
                kind: segment.kind.clone(),
            });
        }

        if segment.skip_mode.is_implicit() {
            for (pos_idx, read_pos) in positions.iter().enumerate() {
                if used_positions.get(pos_idx).copied().unwrap_or(false) {
                    continue;
                }

                calls.push(ModificationCall {
                    position: *read_pos,
                    probability: None,
                    methylated: false,
                    kind: segment.kind.clone(),
                });
            }
        }
    }

    if let Some(values) = ml {
        if probability_index != values.len() {
            bail!(
                "ML tag contains {} values but {} modifications were parsed",
                values.len(),
                probability_index
            );
        }
    }

    Ok(calls)
}

fn canonical_positions(sequence: &[u8], base: u8, reverse: bool) -> Vec<usize> {
    let mut positions: Vec<usize> = sequence
        .iter()
        .enumerate()
        .filter_map(|(idx, b)| (b.to_ascii_uppercase() == base).then_some(idx))
        .collect();

    if reverse {
        positions.reverse();
    }

    positions
}

#[derive(Clone)]
struct Segment {
    canonical_base: u8,
    reverse: bool,
    deltas: Vec<usize>,
    kind: ModificationKind,
    skip_mode: SkipMode,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SkipMode {
    Explicit,
    ImplicitUnmodified,
    DefaultImplicitUnmodified,
}

impl SkipMode {
    fn is_implicit(self) -> bool {
        matches!(
            self,
            SkipMode::ImplicitUnmodified | SkipMode::DefaultImplicitUnmodified
        )
    }
}

fn parse_segments(mm: &str) -> Result<Vec<Segment>> {
    mm.split_terminator(';')
        .filter(|s| !s.is_empty())
        .map(parse_segment)
        .collect()
}

fn parse_segment(segment: &str) -> Result<Segment> {
    let mut parts = segment.split(',');
    let head = parts
        .next()
        .ok_or_else(|| anyhow!("missing MM segment header"))?;

    let deltas = parts
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<usize>()
                .with_context(|| format!("invalid delta `{s}`"))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut chars = head.chars();
    let base = chars
        .next()
        .ok_or_else(|| anyhow!("missing canonical base in MM segment `{segment}`"))?;
    let strand_char = chars
        .next()
        .ok_or_else(|| anyhow!("missing strand in MM segment `{segment}`"))?;
    let mut code: String = chars.collect();
    if code.is_empty() {
        bail!("missing modification code in MM segment `{segment}`");
    }

    let mut skip_mode = SkipMode::DefaultImplicitUnmodified;
    if let Some(last) = code.chars().last() {
        match last {
            '?' => {
                skip_mode = SkipMode::Explicit;
                code.pop();
            }
            '.' => {
                skip_mode = SkipMode::ImplicitUnmodified;
                code.pop();
            }
            _ => {}
        }
    }

    if code.is_empty() {
        bail!("missing modification code in MM segment `{segment}`");
    }

    let strand = match strand_char {
        '+' => Strand::Forward,
        '-' => Strand::Reverse,
        '.' => Strand::Unknown,
        other => bail!("invalid strand `{other}` in MM segment `{segment}`"),
    };

    let label = long_modification_label(base, &code);
    let canonical_base = base.to_ascii_uppercase() as u8;
    let reverse = matches!(strand, Strand::Reverse);

    Ok(Segment {
        canonical_base,
        reverse,
        deltas,
        kind: ModificationKind {
            canonical_base,
            label,
            code,
            strand,
        },
        skip_mode,
    })
}

fn long_modification_label(base: char, code: &str) -> String {
    match (base, code) {
        ('A', "a") => "6mA".to_string(),
        ('C', "m") => "5mC".to_string(),
        ('C', "21839") => "4mC".to_string(),
        ('C', "C") => "4mC".to_string(),
        ('C', "h") => "5hmC".to_string(),
        ('C', "f") => "5fC".to_string(),
        ('C', "g") => "5gmC".to_string(),
        _ => format!("{base}+{code}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_methylated_segment() {
        let calls = parse_modifications("A+a,0;", b"AT", Some(&[128])).unwrap();
        assert_eq!(calls.len(), 1);
        assert!(calls[0].methylated);
        assert_eq!(calls[0].position, 0);
    }

    #[test]
    fn implicit_segments_record_canonical_calls() {
        let calls = parse_modifications("A+a.,0;", b"AAA", Some(&[255])).unwrap();
        assert_eq!(calls.len(), 3);
        let modified = calls.iter().filter(|c| c.methylated).count();
        let canonical = calls.iter().filter(|c| !c.methylated).count();
        assert_eq!((modified, canonical), (1, 2));
    }

    #[test]
    fn explicit_segments_do_not_record_canonical_calls() {
        let calls = parse_modifications("A+a?,0;", b"AAA", Some(&[255])).unwrap();
        assert_eq!(calls.len(), 1);
        assert!(calls[0].methylated);
    }

    #[test]
    fn reverse_segments_record_calls_on_forward_bases() {
        let calls = parse_modifications("A-a,0;", b"GATC", Some(&[200])).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].position, 1);
        assert_eq!(calls[0].kind.strand, Strand::Reverse);
        assert!(calls[0].methylated);
        assert!((calls[0].probability.unwrap() - 200.0 / 255.0).abs() < 1e-6);
    }

    #[test]
    fn mm_segments_error_when_positions_missing() {
        let err = parse_modifications("A+a,1;", b"A", None).unwrap_err();
        assert!(err.to_string().contains("canonical index"));
    }
}
