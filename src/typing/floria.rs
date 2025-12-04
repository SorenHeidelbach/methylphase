use crate::typing::data::{Categories, CategoryConfig, Dataset};
use crate::typing::error::MethylError;
use ndarray::Array2;
use regex::Regex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug)]
struct HapsetBlock {
    _hap_name: String,
    _hap_index: usize,
    snp_start: usize,
    snp_end: usize,
    entries: Vec<FloriaEntry>,
}

#[derive(Debug, Clone)]
struct FloriaEntry {
    id: String,
    start: usize,
    end: usize,
}

/// Result of parsing a Floria haploset file.
pub struct FloriaParsed {
    pub dataset: Dataset,
    pub config: CategoryConfig,
    /// Subcategory labels per category (same ordering as dataset categories).
    pub subcat_names: Vec<Vec<String>>,
}

/// Parse a Floria haploset file and construct Dataset/Config plus subcategory labels.
///
/// Categories are defined by SNP segments built from hapset breakpoints. Each hapset
/// overlapping a segment is a subcategory within that category.
pub fn parse_floria<P: AsRef<Path>>(
    path: P,
) -> Result<FloriaParsed, MethylError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let header_re =
        Regex::new(r"^>HAP(\d+)[^\t]*\t.*SNPRANGE:(\d+)-(\d+)").map_err(|e| {
            MethylError::Parse(format!("invalid regex for header: {}", e))
        })?;

    let mut blocks: Vec<HapsetBlock> = Vec::new();
    let mut current: Option<HapsetBlock> = None;

    for line_res in reader.lines() {
        let line = line_res?;
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if let Some(block) = current.take() {
                blocks.push(block);
            }
            let caps = header_re.captures(&line).ok_or_else(|| {
                MethylError::Parse(format!("failed to parse header line: {}", line))
            })?;
            let hap_index: usize = caps[1].parse().map_err(|e| {
                MethylError::Parse(format!("invalid hap index in header: {}", e))
            })?;
            let snp_start: usize = caps[2].parse().map_err(|e| {
                MethylError::Parse(format!("invalid snp start: {}", e))
            })?;
            let snp_end: usize = caps[3].parse().map_err(|e| {
                MethylError::Parse(format!("invalid snp end: {}", e))
            })?;
            current = Some(HapsetBlock {
                _hap_name: format!("HAP{}", hap_index),
                _hap_index: hap_index,
                snp_start,
                snp_end,
                entries: Vec::new(),
            });
        } else {
            if line.split_whitespace().count() == 1 {
                // Dense allele string (e.g., "010011?1...") covering the full SNPRANGE.
                let alleles = line.trim();
                let block = current
                    .as_mut()
                    .ok_or_else(|| MethylError::Parse("allele line before header".into()))?;
                let expected_len = block.snp_end - block.snp_start + 1;
                if alleles.len() != expected_len {
                    return Err(MethylError::Parse(format!(
                        "allele string length {} does not match SNPRANGE size {}",
                        alleles.len(),
                        expected_len
                    )));
                }
                for (offset, ch) in alleles.chars().enumerate() {
                    if ch == '?' {
                        continue;
                    }
                    let pos = block.snp_start + offset;
                    block.entries.push(FloriaEntry {
                        id: format!("snp_{}", pos),
                        start: pos,
                        end: pos,
                    });
                }
                // If we saw a dense allele string, stop reading further lines until next header.
                // Floria dense format has one line per hapset.
                continue;
            } else {
                let mut parts = line.split_whitespace();
                let id = parts
                    .next()
                    .ok_or_else(|| MethylError::Parse("missing id column".into()))?
                    .to_string();
                let start: usize = parts
                    .next()
                    .ok_or_else(|| MethylError::Parse("missing start column".into()))?
                    .parse()
                    .map_err(|e| MethylError::Parse(format!("invalid start: {}", e)))?;
                let end: usize = parts
                    .next()
                    .ok_or_else(|| MethylError::Parse("missing end column".into()))?
                    .parse()
                    .map_err(|e| MethylError::Parse(format!("invalid end: {}", e)))?;
                if let Some(block) = current.as_mut() {
                    block.entries.push(FloriaEntry { id, start, end });
                } else {
                    return Err(MethylError::Parse(
                        "data line encountered before any header".into(),
                    ));
                }
            }
        }
    }
    if let Some(block) = current {
        blocks.push(block);
    }

    if blocks.is_empty() {
        return Err(MethylError::Parse(
            "no hapset blocks found in floria file".into(),
        ));
    }

    // Build SNP boundaries to create non-overlapping segments where coverage can change.
    let mut boundaries: Vec<usize> = Vec::new();
    for block in &blocks {
        boundaries.push(block.snp_start);
        // end is inclusive; use end + 1 for half-open ranges to build segments.
        boundaries.push(block.snp_end + 1);
    }
    boundaries.sort_unstable();
    boundaries.dedup();

    // Build categories as contiguous SNP segments where at least one hapset overlaps.
    let mut category_defs: Vec<(usize, usize, Vec<&HapsetBlock>)> = Vec::new();
    for window in boundaries.windows(2) {
        let seg_start = window[0];
        let seg_end_exclusive = window[1];
        if seg_start >= seg_end_exclusive {
            continue;
        }
        let seg_end_inclusive = seg_end_exclusive - 1;
        let covering: Vec<&HapsetBlock> = blocks
            .iter()
            .filter(|b| b.snp_start <= seg_start && b.snp_end >= seg_end_inclusive)
            .collect();
        if covering.is_empty() {
            continue;
        }
        category_defs.push((seg_start, seg_end_inclusive, covering));
    }

    // Category names and levels.
    let mut category_names = Vec::new();
    let mut n_levels = Vec::new();
    let mut subcat_names: Vec<Vec<String>> = Vec::new();
    for (start, end, covering) in &category_defs {
        category_names.push(format!("hap_{}_{}", start, end));
        n_levels.push(covering.len());
        subcat_names.push(
            covering
                .iter()
                .map(|b| b._hap_name.clone())
                .collect::<Vec<_>>(),
        );
    }

    // Collect all ids across hapsets.
    let mut id_set = HashMap::<String, usize>::new();
    for block in &blocks {
        for entry in &block.entries {
            let next_idx = id_set.len();
            id_set.entry(entry.id.clone()).or_insert(next_idx);
        }
    }

    let n_rows = id_set.len();
    let n_categories = category_names.len();
    let mut data = Array2::<Option<usize>>::from_elem((n_rows, n_categories), None);

    // Assign entries to categories if their SNP span overlaps the category segment.
    for (cat_idx, (start, end, covering)) in category_defs.iter().enumerate() {
        for (sub_idx, block) in covering.iter().enumerate() {
            for entry in &block.entries {
                if entry.start <= *end && entry.end >= *start {
                    let row = *id_set.get(&entry.id).unwrap();
                    match data[(row, cat_idx)] {
                        Some(existing) if existing != sub_idx => {
                            // Ambiguous assignment; keep the first one.
                        }
                        None => data[(row, cat_idx)] = Some(sub_idx),
                        _ => {}
                    }
                }
            }
        }
    }

    // Stable order ids by insertion.
    let mut ids: Vec<(usize, String)> = id_set
        .into_iter()
        .map(|(id, idx)| (idx, id))
        .collect();
    ids.sort_by_key(|(idx, _)| *idx);
    let ids: Vec<String> = ids.into_iter().map(|(_, id)| id).collect();

    let dataset = Dataset {
        ids,
        data,
        n_levels: n_levels.clone(),
        category_names: category_names.clone(),
        methylation: None,
        methylation_names: Vec::new(),
    };
    dataset.validate()?;

    let config = CategoryConfig {
        categories: Categories {
            names: category_names,
            levels: n_levels,
        },
        methylation_names: Vec::new(),
    };

    Ok(FloriaParsed {
        dataset,
        config,
        subcat_names,
    })
}
