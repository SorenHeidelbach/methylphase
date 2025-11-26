use anyhow::{bail, Context, Result};
use csv::{ReaderBuilder, StringRecord, Trim};
use std::{collections::HashSet, fs, path::Path};

pub fn load_motif_file(path: &Path, allowed_ids: Option<&[String]>) -> Result<Vec<String>> {
    let data =
        fs::read(path).with_context(|| format!("failed to read motif file {}", path.display()))?;
    if data.is_empty() {
        bail!("motif file {} is empty", path.display());
    }

    let delimiter = detect_delimiter(&data);
    let mut reader = ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(true)
        .trim(Trim::All)
        .flexible(true)
        .from_reader(&data[..]);

    let headers = reader
        .headers()
        .with_context(|| format!("failed to read headers from {}", path.display()))?
        .clone();

    let motif_idx = find_column(&headers, "motif")
        .ok_or_else(|| anyhow::anyhow!("motif column missing in {}", path.display()))?;
    let mod_type_idx = find_column(&headers, "mod_type")
        .ok_or_else(|| anyhow::anyhow!("mod_type column missing in {}", path.display()))?;
    let mod_position_idx = find_column(&headers, "mod_position")
        .ok_or_else(|| anyhow::anyhow!("mod_position column missing in {}", path.display()))?;
    let motif_comp_idx = find_column(&headers, "motif_complement");
    let mod_pos_comp_idx = find_column(&headers, "mod_position_complement");
    let bin_idx = ["id", "bin", "reference"]
        .iter()
        .find_map(|name| find_column(&headers, name));

    if motif_comp_idx.is_some() ^ mod_pos_comp_idx.is_some() {
        bail!(
            "motif_complement and mod_position_complement columns must both be present in {}",
            path.display()
        );
    }

    let allowed_lookup = match allowed_ids {
        Some(ids) if !ids.is_empty() => Some(ids.iter().cloned().collect::<HashSet<String>>()),
        _ => None,
    };

    let mut specs = Vec::new();
    for record in reader.records() {
        let record = record
            .with_context(|| format!("failed to read motif file record in {}", path.display()))?;

        if let (Some(idx), Some(ref allowed)) = (bin_idx, allowed_lookup.as_ref()) {
            let identifier = record.get(idx).map(str::trim).unwrap_or("");
            if identifier.is_empty() || !allowed.contains(identifier) {
                continue;
            }
        }

        let motif = record
            .get(motif_idx)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| anyhow::anyhow!("row missing motif value"))?;
        let motif_len = motif.chars().count();

        let mod_label = record
            .get(mod_type_idx)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| anyhow::anyhow!("row missing mod_type value"))
            .and_then(normalize_mod_label)?;

        let position_raw = record
            .get(mod_position_idx)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| anyhow::anyhow!("row missing mod_position value"))?;

        let mod_position: usize = position_raw.parse().with_context(|| {
            format!(
                "invalid mod_position value `{}` in {}",
                position_raw,
                path.display()
            )
        })?;
        if mod_position >= motif_len {
            bail!(
                "mod_position {} exceeds motif length {} in {}",
                mod_position,
                motif_len,
                path.display()
            );
        }

        let spec = format!("{}_{}_{}", motif, mod_label, mod_position);
        specs.push(spec.clone());

        if let (Some(motif_idx), Some(pos_idx)) = (motif_comp_idx, mod_pos_comp_idx) {
            let comp_motif = record
                .get(motif_idx)
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string());
            let comp_position = record.get(pos_idx).map(str::trim).filter(|s| !s.is_empty());

            if let (Some(comp_motif), Some(comp_position)) = (comp_motif, comp_position) {
                let position: usize = comp_position.parse().with_context(|| {
                    format!(
                        "invalid mod_position_complement value `{}` in {}",
                        comp_position,
                        path.display()
                    )
                })?;
                let comp_len = comp_motif.chars().count();
                if position >= comp_len {
                    bail!(
                        "mod_position_complement {} exceeds motif length {} in {}",
                        position,
                        comp_len,
                        path.display()
                    );
                }
                let comp_spec = format!("{}_{}_{}", comp_motif, mod_label, position);
                if comp_spec != spec {
                    specs.push(comp_spec);
                }
            }
        }
    }

    if specs.is_empty() {
        if bin_idx.is_some() && allowed_lookup.is_some() {
            let ids = allowed_lookup
                .as_ref()
                .map(|set| {
                    let mut values: Vec<_> = set.iter().cloned().collect();
                    values.sort();
                    values.join(",")
                })
                .unwrap_or_default();
            bail!(
                "motif file {} produced no entries for ids {}",
                path.display(),
                ids
            );
        } else {
            bail!("motif file {} produced no entries", path.display());
        }
    }

    Ok(specs)
}

fn detect_delimiter(data: &[u8]) -> u8 {
    let first_line = data.split(|b| *b == b'\n').next().unwrap_or(data);
    if first_line.contains(&b'\t') {
        b'\t'
    } else {
        b','
    }
}

fn find_column(headers: &StringRecord, name: &str) -> Option<usize> {
    headers.iter().position(|h| h.eq_ignore_ascii_case(name))
}

fn normalize_mod_label(raw: &str) -> Result<String> {
    let token = raw.trim();
    if token.is_empty() {
        bail!("mod_type value is empty");
    }

    if token.len() == 1 {
        if let Some(label) = code_to_label(token.chars().next().unwrap()) {
            return Ok(label.to_string());
        }
    }

    let lowered = token.to_ascii_lowercase();
    if let Some(label) = name_to_label(&lowered) {
        return Ok(label.to_string());
    }

    Ok(token.to_string())
}

fn code_to_label(code: char) -> Option<&'static str> {
    match code.to_ascii_lowercase() {
        'a' => Some("6mA"),
        'm' => Some("5mC"),
        'h' => Some("5hmC"),
        'f' => Some("5fC"),
        'g' => Some("5gmC"),
        'c' => Some("4mC"),
        _ => None,
    }
}

fn name_to_label(name: &str) -> Option<&'static str> {
    match name {
        "6ma" => Some("6mA"),
        "5mc" => Some("5mC"),
        "5hmc" => Some("5hmC"),
        "5fc" => Some("5fC"),
        "5gmc" => Some("5gmC"),
        "4mc" => Some("4mC"),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_load_motif_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "motif\tmod_type\tmod_position\nGATC\t6mA\t1\nCCWGG\tm\t0"
        )
        .unwrap();
        let specs = load_motif_file(file.path(), None).unwrap();
        assert_eq!(
            specs,
            vec!["GATC_6mA_1".to_string(), "CCWGG_5mC_0".to_string()]
        );
    }

    #[test]
    fn test_load_motif_file_with_complements() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "motif\tmod_type\tmod_position\tmotif_complement\tmod_position_complement\n\
             GATC\t6mA\t1\tGATC\t2\n\
             CCWGG\tm\t0\tCCWGG\t4"
        )
        .unwrap();
        let specs = load_motif_file(file.path(), None).unwrap();
        assert_eq!(
            specs,
            vec![
                "GATC_6mA_1".to_string(),
                "GATC_6mA_2".to_string(),
                "CCWGG_5mC_0".to_string(),
                "CCWGG_5mC_4".to_string(),
            ]
        );
    }

    #[test]
    fn test_load_motif_file_filters_by_bin_column() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "bin\tmotif\tmod_type\tmod_position\n\
             binA\tGATC\t6mA\t1\n\
             binB\tCCWGG\tm\t0"
        )
        .unwrap();
        let ids = vec!["binB".to_string()];
        let specs = load_motif_file(file.path(), Some(&ids)).unwrap();
        assert_eq!(specs, vec!["CCWGG_5mC_0".to_string()]);
    }

    #[test]
    fn test_load_motif_file_supports_reference_column_case_insensitive() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "Reference\tmotif\tmod_type\tmod_position\n\
             BINA\tGATC\t6mA\t1\n\
             BINB\tCCWGG\tm\t0"
        )
        .unwrap();
        let ids = vec!["BINB".to_string()];
        let specs = load_motif_file(file.path(), Some(&ids)).unwrap();
        assert_eq!(specs, vec!["CCWGG_5mC_0".to_string()]);
    }

    #[test]
    fn test_load_motif_file_errors_when_bin_missing() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "bin\tmotif\tmod_type\tmod_position\n\
             binA\tGATC\t6mA\t1"
        )
        .unwrap();
        let ids = vec!["binB".to_string()];
        let err = load_motif_file(file.path(), Some(&ids)).unwrap_err();
        assert!(
            err.to_string().contains("produced no entries for ids binB"),
            "unexpected error: {err}"
        );
    }
}
