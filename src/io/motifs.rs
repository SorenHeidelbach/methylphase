use anyhow::{bail, Context, Result};
use csv::{ReaderBuilder, StringRecord, Trim};
use std::{fs, path::Path};

pub fn load_motif_file(path: &Path) -> Result<Vec<String>> {
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

    let mut specs = Vec::new();
    for record in reader.records() {
        let record = record
            .with_context(|| format!("failed to read motif file record in {}", path.display()))?;

        let motif = record
            .get(motif_idx)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| anyhow::anyhow!("row missing motif value"))?;

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

        specs.push(format!("{}_{}_{}", motif, mod_label, mod_position));
    }

    if specs.is_empty() {
        bail!("motif file {} produced no entries", path.display());
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
            "motif\tmod_type\tmod_position\nGATC\t6mA\t2\nCCWGG\tm\t1"
        )
        .unwrap();
        let specs = load_motif_file(file.path()).unwrap();
        assert_eq!(
            specs,
            vec!["GATC_6mA_2".to_string(), "CCWGG_5mC_1".to_string()]
        );
    }
}
