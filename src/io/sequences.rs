use anyhow::{anyhow, Context, Result};
use bstr::ByteSlice;
use noodles_bam as bam;
use noodles_fastq as fastq;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Seek, SeekFrom},
    path::{Path, PathBuf},
    sync::Mutex,
};

/// Lookup table for read sequences sourced from auxiliary files.
pub enum SequenceCache {
    Memory(HashMap<String, Vec<u8>>),
    FastqIndexed(FastqIndexed),
}

impl SequenceCache {
    pub fn from_path<P: AsRef<Path>>(path: P, index: Option<&Path>) -> Result<Self> {
        let path = path.as_ref();
        match detect_format(path)? {
            SequenceFormat::Bam => Self::from_bam(path),
            SequenceFormat::Fastq => {
                if let Some(idx) = index {
                    Self::from_fastq_indexed(path, idx)
                } else {
                    Self::from_fastq(path)
                }
            }
        }
    }

    pub fn sequence(&self, read_id: &str) -> Option<Vec<u8>> {
        match self {
            SequenceCache::Memory(map) => map.get(read_id).cloned(),
            SequenceCache::FastqIndexed(indexed) => indexed.sequence(read_id),
        }
    }

    fn from_bam(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("failed to open BAM {}", path.display()))?;
        let mut reader = bam::io::Reader::new(file);
        reader
            .read_header()
            .context("failed to read auxiliary BAM header")?;
        let mut record = bam::Record::default();
        let mut sequences = HashMap::new();
        let mut count = 0usize;
        eprintln!(
            "sequence-fallback: loading sequences from BAM {}...",
            path.display()
        );
        loop {
            if reader.read_record(&mut record)? == 0 {
                break;
            }
            let name = record
                .name()
                .context("encountered unnamed record in auxiliary BAM")?
                .to_str()
                .context("auxiliary BAM read name contains invalid UTF-8")?
                .to_string();
            let seq: Vec<u8> = record.sequence().iter().collect();
            sequences.insert(name, seq);
            count += 1;
            if count % 100_000 == 0 {
                eprintln!(
                    "sequence-fallback: loaded {} sequences from {}...",
                    count,
                    path.display()
                );
            }
        }
        eprintln!(
            "sequence-fallback: finished loading {} sequences from {}",
            count,
            path.display()
        );
        Ok(Self::Memory(sequences))
    }

    fn from_fastq(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("failed to open FASTQ {}", path.display()))?;
        let mut reader = fastq::io::Reader::new(BufReader::new(file));
        let mut record = fastq::Record::default();
        let mut sequences = HashMap::new();
        let mut count = 0usize;
        eprintln!(
            "sequence-fallback: loading sequences from FASTQ {}...",
            path.display()
        );
        loop {
            if reader.read_record(&mut record)? == 0 {
                break;
            }
            let definition = record.definition();
            let name = definition
                .name()
                .to_str()
                .context("FASTQ read name contains invalid UTF-8")?
                .to_string();
            if name.is_empty() {
                continue;
            }
            let seq = record.sequence().to_vec();
            sequences.insert(name, seq);
            count += 1;
            if count % 100_000 == 0 {
                eprintln!(
                    "sequence-fallback: loaded {} sequences from {}...",
                    count,
                    path.display()
                );
            }
        }
        eprintln!(
            "sequence-fallback: finished loading {} sequences from {}",
            count,
            path.display()
        );
        Ok(Self::Memory(sequences))
    }

    fn from_fastq_indexed(path: &Path, index_path: &Path) -> Result<Self> {
        eprintln!(
            "sequence-fallback: using FASTQ index {} for {}",
            index_path.display(),
            path.display()
        );
        let entries = load_fastq_index(index_path)?;
        let reader =
            File::open(path).with_context(|| format!("failed to open FASTQ {}", path.display()))?;
        Ok(Self::FastqIndexed(FastqIndexed {
            reader: Mutex::new(BufReader::new(reader)),
            entries,
            path: path.to_path_buf(),
        }))
    }
}

pub struct FastqIndexed {
    reader: Mutex<BufReader<File>>,
    entries: HashMap<String, u64>,
    path: PathBuf,
}

impl FastqIndexed {
    fn sequence(&self, read_id: &str) -> Option<Vec<u8>> {
        let offset = *self.entries.get(read_id)?;
        let mut reader = self.reader.lock().ok()?;
        reader.seek(SeekFrom::Start(offset)).ok()?;
        let mut header = String::new();
        reader.read_line(&mut header).ok()?;
        if !header.starts_with('@') {
            eprintln!(
                "sequence-fallback: malformed FASTQ record for {} in {}",
                read_id,
                self.path.display()
            );
            return None;
        }
        let mut seq_line = String::new();
        reader.read_line(&mut seq_line).ok()?;
        let sequence = seq_line
            .trim_end_matches(|c| c == '\r' || c == '\n')
            .as_bytes()
            .to_vec();
        let mut plus = String::new();
        reader.read_line(&mut plus).ok()?;
        let mut qual = String::new();
        reader.read_line(&mut qual).ok()?;
        Some(sequence)
    }
}

fn load_fastq_index(path: &Path) -> Result<HashMap<String, u64>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open FASTQ index {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut entries = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut parts = trimmed.split('\t');
        let read_id = parts
            .next()
            .ok_or_else(|| anyhow!("invalid fqidx line `{line}`"))?;
        let offset_str = parts
            .next()
            .ok_or_else(|| anyhow!("invalid fqidx line `{line}`"))?;
        let offset = offset_str
            .parse::<u64>()
            .with_context(|| format!("invalid offset `{offset_str}` in fqidx line `{line}`"))?;
        entries.insert(read_id.to_string(), offset);
    }
    eprintln!(
        "sequence-fallback: loaded {} index entries from {}",
        entries.len(),
        path.display()
    );
    Ok(entries)
}

enum SequenceFormat {
    Bam,
    Fastq,
}

fn detect_format(path: &Path) -> Result<SequenceFormat> {
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext.to_ascii_lowercase().as_str() {
            "bam" => return Ok(SequenceFormat::Bam),
            "fastq" | "fq" => return Ok(SequenceFormat::Fastq),
            _ => {}
        }
    }

    let filename = path
        .file_name()
        .and_then(|f| f.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if filename.ends_with(".fastq") || filename.ends_with(".fq") {
        Ok(SequenceFormat::Fastq)
    } else if filename.ends_with(".bam") {
        Ok(SequenceFormat::Bam)
    } else {
        Err(anyhow!(
            "unable to infer sequence file format for {} (expected .bam or .fastq)",
            path.display()
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::Builder;

    #[test]
    fn loads_sequences_from_fastq() -> Result<()> {
        let mut tmp = Builder::new().suffix(".fastq").tempfile()?;
        std::io::Write::write_all(&mut tmp, b"@read0\nACGT\n+\n!!!!\n@read1\nTTAA\n+\n####\n")?;
        let cache = SequenceCache::from_path(tmp.path(), None)?;
        assert_eq!(cache.sequence("read0"), Some(b"ACGT".to_vec()));
        assert_eq!(cache.sequence("read1"), Some(b"TTAA".to_vec()));
        assert!(cache.sequence("nope").is_none());
        Ok(())
    }
}
