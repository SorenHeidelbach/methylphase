use anyhow::{bail, Context, Result};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

#[derive(Debug, Clone)]
pub struct ContigBins {
    bins: Vec<BinAssignment>,
}

#[derive(Debug, Clone)]
pub struct BinAssignment {
    name: String,
    safe_name: String,
    contigs: Vec<String>,
}

impl ContigBins {
    pub fn from_path(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("failed to open contig bins {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut bins = Vec::new();
        let mut index = HashMap::new();
        let mut safe_counts: HashMap<String, usize> = HashMap::new();

        for (line_idx, line) in reader.lines().enumerate() {
            let line = line.with_context(|| {
                format!(
                    "failed to read contig bins line {} from {}",
                    line_idx + 1,
                    path.display()
                )
            })?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            let mut splitter = trimmed.splitn(2, char::is_whitespace);
            let contig = match splitter.next() {
                Some(value) if !value.is_empty() => value,
                _ => continue,
            };
            let bin = match splitter.next() {
                Some(value) if !value.trim().is_empty() => value.trim(),
                _ => {
                    bail!(
                        "contig bins file {} line {} is missing the bin column",
                        path.display(),
                        line_idx + 1
                    );
                }
            };

            if contig.eq_ignore_ascii_case("contig") && bin.eq_ignore_ascii_case("bin") {
                continue;
            }

            let contig_name = contig.to_string();
            let bin_name = bin.to_string();
            let entry_idx = *index.entry(bin_name.clone()).or_insert_with(|| {
                let safe = unique_safe_name(&bin_name, &mut safe_counts);
                bins.push(BinAssignment::new(bin_name.clone(), safe));
                bins.len() - 1
            });
            bins[entry_idx].contigs.push(contig_name);
        }

        if bins.is_empty() {
            bail!(
                "no contig bins were found in {}; expected '<contig>\\t<bin>' rows",
                path.display()
            );
        }

        Ok(Self { bins })
    }

    pub fn iter(&self) -> impl Iterator<Item = &BinAssignment> {
        self.bins.iter()
    }
}

impl BinAssignment {
    fn new(name: String, safe_name: String) -> Self {
        Self {
            name,
            safe_name,
            contigs: Vec::new(),
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn safe_name(&self) -> &str {
        &self.safe_name
    }

    pub fn contigs(&self) -> &[String] {
        &self.contigs
    }
}

fn unique_safe_name(name: &str, counts: &mut HashMap<String, usize>) -> String {
    let mut safe = sanitize_name(name);
    let counter = counts.entry(safe.clone()).or_insert(0);
    if *counter > 0 {
        safe.push('_');
        safe.push_str(&counter.to_string());
    }
    *counter += 1;
    safe
}

fn sanitize_name(name: &str) -> String {
    let mut sanitized = String::with_capacity(name.len());
    for ch in name.chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
            sanitized.push(ch);
        } else {
            sanitized.push('_');
        }
    }
    if sanitized.is_empty() {
        "_".to_string()
    } else {
        sanitized
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn parses_bins_and_sanitizes_names() -> Result<()> {
        let mut tmp = NamedTempFile::new()?;
        writeln!(tmp, "# comment")?;
        writeln!(tmp, "contig\tbin")?;
        writeln!(tmp, "chr1\tbin 1")?;
        writeln!(tmp, "chr2\tbin 1")?;
        writeln!(tmp, "chr3\tweird/bin")?;
        writeln!(tmp, "chr4\tweird/bin")?;

        let bins = ContigBins::from_path(tmp.path())?;
        let assignments: Vec<_> = bins.iter().collect();
        assert!(!assignments.is_empty());
        assert_eq!(assignments.len(), 2);
        assert_eq!(assignments[0].name(), "bin 1");
        assert_eq!(assignments[0].safe_name(), "bin_1");
        assert_eq!(
            assignments[0].contigs(),
            &["chr1".to_string(), "chr2".to_string()]
        );
        assert_eq!(assignments[1].name(), "weird/bin");
        assert!(assignments[1].safe_name().starts_with("weird_bin"));
        assert_eq!(assignments[1].contigs().len(), 2);
        Ok(())
    }
}
