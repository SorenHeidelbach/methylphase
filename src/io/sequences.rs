use anyhow::{anyhow, bail, Context, Result};
use bstr::ByteSlice;
use flate2::read::MultiGzDecoder;
use noodles_bam as bam;
use noodles_bgzf::{self as bgzf, io::Seek as _};
use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    sync::Mutex,
};
use tempfile::NamedTempFile;

/// Lookup table for read sequences sourced from auxiliary files.
pub enum SequenceCache {
    FastqIndexed(FastqIndexed),
    BamIndexed(BamIndexed),
}

impl SequenceCache {
    pub fn from_path<P: AsRef<Path>>(path: P, index: Option<&Path>) -> Result<Self> {
        let path = path.as_ref();
        match detect_format(path)? {
            SequenceFormat::Bam => Self::from_bam_with_index(path, index),
            SequenceFormat::Fastq => Self::from_fastq_with_index(path, index),
            SequenceFormat::FastqGz => Self::from_gz_fastq(path, index),
        }
    }

    pub fn sequence(&self, read_id: &str) -> Option<Vec<u8>> {
        match self {
            SequenceCache::FastqIndexed(indexed) => indexed.sequence(read_id),
            SequenceCache::BamIndexed(indexed) => indexed.sequence(read_id),
        }
    }

    fn from_bam_with_index(path: &Path, index: Option<&Path>) -> Result<Self> {
        let index_path = index
            .map(PathBuf::from)
            .unwrap_or_else(|| default_index_name(path, "bai"));

        let entries = if index_path.exists() {
            load_bam_index(&index_path)?
        } else {
            build_bam_index(path, &index_path)?
        };

        let file =
            File::open(path).with_context(|| format!("failed to open BAM {}", path.display()))?;
        let mut reader = bam::io::Reader::new(file);
        reader
            .read_header()
            .context("failed to read BAM header for sequence fallback")?;
        Ok(Self::BamIndexed(BamIndexed {
            reader: Mutex::new(reader),
            entries,
        }))
    }

    fn from_fastq_with_index(path: &Path, index_hint: Option<&Path>) -> Result<Self> {
        let fastq_path = path.to_path_buf();
        let (index_path, entries) = resolve_fastq_index(&fastq_path, index_hint)?;
        eprintln!(
            "sequence-fallback: using FASTQ index {} for {} ({} entries)",
            index_path.display(),
            fastq_path.display(),
            entries.len()
        );
        let reader = File::open(&fastq_path)
            .with_context(|| format!("failed to open FASTQ {}", fastq_path.display()))?;
        Ok(Self::FastqIndexed(FastqIndexed {
            reader: Mutex::new(BufReader::new(reader)),
            entries,
            path: fastq_path,
        }))
    }

    fn from_gz_fastq(path: &Path, index_hint: Option<&Path>) -> Result<Self> {
        let materialized = materialize_fastq_from_gzip(path)?;
        Self::from_fastq_with_index(&materialized, index_hint)
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

pub struct BamIndexed {
    reader: Mutex<bam::io::Reader<bgzf::io::Reader<File>>>,
    entries: HashMap<String, bgzf::VirtualPosition>,
}

impl BamIndexed {
    fn sequence(&self, read_id: &str) -> Option<Vec<u8>> {
        let position = *self.entries.get(read_id)?;
        let mut reader = self.reader.lock().ok()?;
        reader.get_mut().seek_to_virtual_position(position).ok()?;
        let mut record = bam::Record::default();
        if reader.read_record(&mut record).ok()? == 0 {
            return None;
        }
        let seq: Vec<u8> = record.sequence().iter().collect();
        Some(seq)
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
    Ok(entries)
}

fn build_fastq_index(fastq_path: &Path, index_path: &Path) -> Result<HashMap<String, u64>> {
    eprintln!(
        "sequence-fallback: building FASTQ index {} for {}",
        index_path.display(),
        fastq_path.display()
    );
    let file = File::open(fastq_path)
        .with_context(|| format!("failed to open FASTQ {}", fastq_path.display()))?;
    let mut reader = BufReader::new(file);
    let mut entries = HashMap::new();
    let mut offset: u64 = 0;
    loop {
        let record_offset = offset;
        let mut header = String::new();
        let bytes = reader.read_line(&mut header)?;
        if bytes == 0 {
            break;
        }
        offset += bytes as u64;
        if !header.starts_with('@') {
            bail!(
                "invalid FASTQ header `{}` while building index for {}",
                header.trim_end(),
                fastq_path.display()
            );
        }
        let name = header
            .trim_end_matches(|c| c == '\r' || c == '\n')
            .get(1..)
            .and_then(|rest| rest.split_whitespace().next())
            .ok_or_else(|| anyhow!("invalid FASTQ header `{header}`"))?
            .to_string();
        if name.is_empty() {
            bail!(
                "encountered FASTQ record without a name in {}",
                fastq_path.display()
            );
        }

        let mut seq = String::new();
        offset += reader.read_line(&mut seq)? as u64;
        let mut plus = String::new();
        let plus_bytes = reader.read_line(&mut plus)?;
        if plus_bytes == 0 || !plus.starts_with('+') {
            bail!(
                "invalid FASTQ separator for read {} in {}",
                name,
                fastq_path.display()
            );
        }
        offset += plus_bytes as u64;
        let mut qual = String::new();
        let qual_bytes = reader.read_line(&mut qual)?;
        if qual_bytes == 0 {
            bail!(
                "missing FASTQ quality line for read {} in {}",
                name,
                fastq_path.display()
            );
        }
        offset += qual_bytes as u64;
        entries.insert(name, record_offset);
    }

    let tmp_ext = index_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|ext| format!("{ext}.tmp"))
        .unwrap_or_else(|| "tmp".to_string());
    let mut tmp_path = index_path.to_path_buf();
    tmp_path.set_extension(tmp_ext);
    let mut writer = BufWriter::new(
        File::create(&tmp_path)
            .with_context(|| format!("failed to create {}", tmp_path.display()))?,
    );
    let mut sorted: Vec<_> = entries.iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(b.0));
    for (read_id, offset) in sorted {
        writeln!(writer, "{read_id}\t{offset}")?;
    }
    writer.flush()?;
    fs::rename(&tmp_path, index_path)
        .with_context(|| format!("failed to write FASTQ index {}", index_path.display()))?;
    eprintln!(
        "sequence-fallback: finished building FASTQ index {} ({} entries)",
        index_path.display(),
        entries.len()
    );
    Ok(entries)
}

fn resolve_fastq_index(
    path: &Path,
    index_hint: Option<&Path>,
) -> Result<(PathBuf, HashMap<String, u64>)> {
    if let Some(explicit) = index_hint {
        if explicit.exists() {
            let entries = load_fastq_index(explicit)?;
            return Ok((explicit.to_path_buf(), entries));
        }
        let entries = build_fastq_index(path, explicit)?;
        return Ok((explicit.to_path_buf(), entries));
    }

    let default = default_index_name(path, "fqidx");
    if default.exists() {
        let entries = load_fastq_index(&default)?;
        Ok((default, entries))
    } else {
        let entries = build_fastq_index(path, &default)?;
        Ok((default, entries))
    }
}

fn default_index_name(path: &Path, ext: &str) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(format!(".{ext}"));
    PathBuf::from(name)
}

fn materialize_fastq_from_gzip(path: &Path) -> Result<PathBuf> {
    let mut target = PathBuf::from(path);
    target.set_extension("");
    if target == path {
        target = path.with_extension("fastq");
    }
    if target.exists() {
        return Ok(target);
    }

    eprintln!(
        "sequence-fallback: decompressing {} to {}",
        path.display(),
        target.display()
    );
    let source = File::open(path)
        .with_context(|| format!("failed to open gzipped FASTQ {}", path.display()))?;
    let mut decoder = MultiGzDecoder::new(source);
    let tmp = NamedTempFile::new_in(
        target
            .parent()
            .ok_or_else(|| anyhow!("unable to determine directory for {}", target.display()))?,
    )
    .context("failed to create temporary FASTQ during decompression")?;
    let mut writer = tmp
        .reopen()
        .context("failed to reopen temporary FASTQ during decompression")?;
    std::io::copy(&mut decoder, &mut writer)
        .context("failed to decompress gzipped FASTQ fallback")?;
    let temp_path = tmp.into_temp_path();
    temp_path.persist(&target).map_err(|e| {
        anyhow!(
            "failed to persist decompressed FASTQ {}: {}",
            target.display(),
            e
        )
    })?;
    Ok(target)
}

fn load_bam_index(path: &Path) -> Result<HashMap<String, bgzf::VirtualPosition>> {
    let file =
        File::open(path).with_context(|| format!("failed to open BAM index {}", path.display()))?;
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
            .ok_or_else(|| anyhow!("invalid BAM index line `{line}`"))?;
        let offset_str = parts
            .next()
            .ok_or_else(|| anyhow!("invalid BAM index line `{line}`"))?;
        let offset = offset_str
            .parse::<u64>()
            .with_context(|| format!("invalid offset `{offset_str}` in BAM index line `{line}`"))?;
        entries.insert(read_id.to_string(), bgzf::VirtualPosition::from(offset));
    }
    Ok(entries)
}

fn build_bam_index(
    bam_path: &Path,
    index_path: &Path,
) -> Result<HashMap<String, bgzf::VirtualPosition>> {
    eprintln!(
        "sequence-fallback: building BAM index {} for {}",
        index_path.display(),
        bam_path.display()
    );
    let file = File::open(bam_path)
        .with_context(|| format!("failed to open BAM {}", bam_path.display()))?;
    let mut reader = bam::io::Reader::new(file);
    reader
        .read_header()
        .context("failed to read BAM header while building sequence index")?;
    let mut record = bam::Record::default();
    let mut entries = HashMap::new();
    loop {
        let pos = reader.get_ref().virtual_position();
        if reader.read_record(&mut record)? == 0 {
            break;
        }
        let Some(name) = record.name().map(|n| n.to_owned()) else {
            continue;
        };
        let read_id = name
            .to_str()
            .context("BAM sequence fallback read name contains invalid UTF-8")?
            .to_string();
        entries.entry(read_id).or_insert(pos);
    }

    let mut writer = BufWriter::new(
        File::create(index_path)
            .with_context(|| format!("failed to create BAM index {}", index_path.display()))?,
    );
    let mut rows: Vec<_> = entries.iter().collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));
    for (read_id, offset) in rows {
        writeln!(writer, "{read_id}\t{}", u64::from(*offset))?;
    }
    writer.flush()?;
    Ok(entries)
}

enum SequenceFormat {
    Bam,
    Fastq,
    FastqGz,
}

fn detect_format(path: &Path) -> Result<SequenceFormat> {
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext.to_ascii_lowercase().as_str() {
            "bam" => return Ok(SequenceFormat::Bam),
            "fastq" | "fq" => return Ok(SequenceFormat::Fastq),
            "gz" => {
                if let Some(stem_ext) = path
                    .file_stem()
                    .and_then(|s| Path::new(s).extension())
                    .and_then(|s| s.to_str())
                {
                    let lower = stem_ext.to_ascii_lowercase();
                    if lower == "fastq" || lower == "fq" {
                        return Ok(SequenceFormat::FastqGz);
                    }
                }
            }
            _ => {}
        }
    }

    let filename = path
        .file_name()
        .and_then(|f| f.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if filename.ends_with(".fastq.gz") || filename.ends_with(".fq.gz") {
        Ok(SequenceFormat::FastqGz)
    } else if filename.ends_with(".fastq") || filename.ends_with(".fq") {
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
    use flate2::{write::GzEncoder, Compression};
    use noodles_bam as bam;
    use noodles_sam::{self as sam, alignment::io::Write as SamAlignmentWrite};
    use std::{
        fs,
        io::{Cursor, Write},
    };
    use tempfile::Builder;

    #[test]
    fn loads_sequences_from_fastq() -> Result<()> {
        let mut tmp = Builder::new().suffix(".fastq").tempfile()?;
        Write::write_all(&mut tmp, b"@read0\nACGT\n+\n!!!!\n@read1\nTTAA\n+\n####\n")?;
        let cache = SequenceCache::from_path(tmp.path(), None)?;
        assert_eq!(cache.sequence("read0"), Some(b"ACGT".to_vec()));
        assert_eq!(cache.sequence("read1"), Some(b"TTAA".to_vec()));
        assert!(cache.sequence("nope").is_none());
        let mut index_name = tmp.path().as_os_str().to_os_string();
        index_name.push(".fqidx");
        let index_path = PathBuf::from(index_name);
        assert!(index_path.exists());
        Ok(())
    }

    #[test]
    fn builds_fastq_index_when_missing() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let fastq_path = dir.path().join("reads.fastq");
        fs::write(
            &fastq_path,
            "@r0\nACGT\n+\n!!!!\n@r1 extra\nTTAA\n+\n####\n",
        )?;
        let cache = SequenceCache::from_path(&fastq_path, None)?;
        assert_eq!(cache.sequence("r0"), Some(b"ACGT".to_vec()));
        assert_eq!(cache.sequence("r1"), Some(b"TTAA".to_vec()));
        let mut index_name = fastq_path.as_os_str().to_os_string();
        index_name.push(".fqidx");
        let index_path = PathBuf::from(index_name);
        assert!(index_path.exists());
        Ok(())
    }

    #[test]
    fn loads_sequences_from_gz_fastq() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let gz_path = dir.path().join("reads.fastq.gz");
        {
            let file = File::create(&gz_path)?;
            let mut encoder = GzEncoder::new(file, Compression::default());
            encoder.write_all(b"@r0\nACGT\n+\n!!!!\n@r1\nTTAA\n+\n####\n")?;
            encoder.finish()?;
        }
        let cache = SequenceCache::from_path(&gz_path, None)?;
        assert_eq!(cache.sequence("r0"), Some(b"ACGT".to_vec()));
        assert_eq!(cache.sequence("r1"), Some(b"TTAA".to_vec()));
        let mut materialized = gz_path.clone();
        materialized.set_extension("");
        if materialized == gz_path {
            materialized = gz_path.with_extension("fastq");
        }
        let mut index_name = materialized.as_os_str().to_os_string();
        index_name.push(".fqidx");
        assert!(PathBuf::from(index_name).exists());
        Ok(())
    }

    #[test]
    fn loads_sequences_from_bam_with_index() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bam_path = dir.path().join("reads.bam");
        write_test_bam(&bam_path)?;

        let cache = SequenceCache::from_path(&bam_path, None)?;
        assert_eq!(cache.sequence("read0"), Some(b"ACGT".to_vec()));
        let mut index_name = bam_path.as_os_str().to_os_string();
        index_name.push(".bai");
        assert!(PathBuf::from(index_name).exists());
        Ok(())
    }

    fn write_test_bam(path: &Path) -> Result<()> {
        const SAM_TEMPLATE: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ctg\tLN:50
read0\t0\tctg\t1\t60\t4M\t*\t0\t0\tACGT\tFFFF
";

        let mut reader = sam::io::Reader::new(Cursor::new(SAM_TEMPLATE.as_bytes()));
        let header = reader.read_header()?;
        let mut record = sam::Record::default();

        let file = File::create(path)?;
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header)?;

        while reader.read_record(&mut record)? != 0 {
            writer.write_alignment_record(&header, &record)?;
        }

        writer.try_finish()?;
        Ok(())
    }
}
