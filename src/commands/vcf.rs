use crate::{
    commands::shared,
    core::{ModificationCall, MotifQuery, ReadRecord},
    io::{BamReader, BamRecordDecoder, SequenceCache},
};
use anyhow::{bail, Context, Result};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    process::Command,
    sync::Arc,
};

const FORMAT_FIELDS: &str = "DP:METH:UNMETH:FRAC:MEANP";

pub fn run(
    bam: PathBuf,
    motifs: Vec<String>,
    motif_file: Option<PathBuf>,
    sequence_fallback: Option<PathBuf>,
    sequence_index: Option<PathBuf>,
    output: PathBuf,
    sample_name: Option<String>,
    methylation_threshold: f32,
    plain_output: bool,
    bin_ids: Option<Vec<String>>,
    contigs: Vec<String>,
) -> Result<()> {
    let motif_queries = shared::load_motif_queries(motifs, motif_file, bin_ids.as_deref())?;
    let sample_name = sample_name.unwrap_or_else(|| "sample".to_string());
    let output_paths = resolve_output_paths(output, !plain_output);

    let mut reader =
        BamReader::open(&bam).with_context(|| format!("failed to open BAM {}", bam.display()))?;
    let sequence_cache = sequence_fallback
        .as_ref()
        .map(|path| {
            SequenceCache::from_path(path, sequence_index.as_deref())
                .with_context(|| format!("failed to load sequence fallback {}", path.display()))
        })
        .transpose()?
        .map(Arc::new);
    let decoder = BamRecordDecoder::new(reader.header(), sequence_cache);

    let contig_filter = if contigs.is_empty() {
        None
    } else {
        Some(contigs)
    };
    eprintln!(
        "vcf: processing BAM {} with {} motifs",
        bam.display(),
        motif_queries.len()
    );

    let mut table = SiteTable::new();
    let mut processed_reads = 0usize;
    reader.visit_records(contig_filter.as_ref().map(|c| c.as_slice()), |record| {
        let read = decoder.decode(record)?;
        process_read(&read, &motif_queries, methylation_threshold, &mut table);
        processed_reads += 1;
        if processed_reads % 1_000 == 0 {
            eprintln!("vcf: processed {} reads...", processed_reads);
        }
        Ok(())
    })?;

    decoder.report();

    if let Some(parent) = output_paths.write_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
    }

    write_vcf(
        &output_paths.write_path,
        &table,
        &motif_queries,
        &sample_name,
    )?;

    if output_paths.compress {
        let gz_path = bgzip_file(&output_paths.write_path)?;
        tabix_vcf(&gz_path)?;
        eprintln!(
            "vcf: wrote {} records to {} (processed {} reads)",
            table.len(),
            gz_path.display(),
            processed_reads
        );
    } else {
        eprintln!(
            "vcf: wrote {} records to {} (processed {} reads)",
            table.len(),
            output_paths.write_path.display(),
            processed_reads
        );
    }

    Ok(())
}

fn process_read(read: &ReadRecord, motifs: &[MotifQuery], threshold: f32, table: &mut SiteTable) {
    let Some(contig) = read.contig_name() else {
        return;
    };

    for (idx, motif) in motifs.iter().enumerate() {
        let motif_offset = motif.mod_position();
        for start in motif.matches(read.sequence()) {
            let mod_pos = start + motif_offset;
            let Some(reference_position) = read.reference_position(mod_pos) else {
                continue;
            };
            let Some(call) = read.modification_at(mod_pos, motif.mod_label()) else {
                continue;
            };
            let methylated = classify_call(call, threshold);
            table.record_call(
                contig,
                reference_position,
                idx,
                methylated,
                call.probability,
            );
        }
    }
}

fn classify_call(call: &ModificationCall, threshold: f32) -> bool {
    if !call.methylated {
        return false;
    }

    if let Some(prob) = call.probability {
        prob >= threshold
    } else {
        true
    }
}

fn write_vcf(
    output: &PathBuf,
    table: &SiteTable,
    motifs: &[MotifQuery],
    sample_name: &str,
) -> Result<()> {
    let file =
        File::create(output).with_context(|| format!("unable to create {}", output.display()))?;
    let mut writer = BufWriter::new(file);
    write_headers(&mut writer, table, sample_name)?;

    for (key, stats) in table.iter() {
        let motif = &motifs[key.motif_idx];
        let depth = stats.depth();
        if depth == 0 {
            continue;
        }
        let ref_idx = motif.mod_position();
        let ref_base = motif
            .motif()
            .get(ref_idx)
            .copied()
            .map(char::from)
            .unwrap_or('N');
        let alt = format_alt(motif.mod_label());
        let fraction = stats.methylation_fraction();
        let info = format!(
            "MOTIF={};MODL={};DP={};METH={};UNMETH={};FRAC={:.4}",
            motif.raw(),
            motif.mod_label(),
            depth,
            stats.methylated,
            stats.unmethylated,
            fraction
        );
        let mean_prob = stats
            .mean_probability()
            .map(|v| format!("{:.4}", v))
            .unwrap_or_else(|| ".".to_string());
        let format_values = format!(
            "{}:{}:{}:{:.4}:{}",
            depth, stats.methylated, stats.unmethylated, fraction, mean_prob
        );

        writeln!(
            writer,
            "{}\t{}\t.\t{}\t{}\t.\tPASS\t{}\t{}\t{}",
            key.contig, key.position, ref_base, alt, info, FORMAT_FIELDS, format_values
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn write_headers<W: Write>(writer: &mut W, table: &SiteTable, sample_name: &str) -> Result<()> {
    writeln!(writer, "##fileformat=VCFv4.3")?;
    writeln!(writer, "##source=methylphase")?;
    for contig in table.contigs() {
        writeln!(writer, "##contig=<ID={contig}>")?;
    }
    writeln!(
        writer,
        "##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"Motif descriptor used for aggregation\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=MODL,Number=1,Type=String,Description=\"Modification label\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total reads covering the motif site\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=METH,Number=1,Type=Integer,Description=\"Reads classified as methylated\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=UNMETH,Number=1,Type=Integer,Description=\"Reads classified as unmethylated\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=FRAC,Number=1,Type=Float,Description=\"Methylated read fraction\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total reads covering the motif site\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=METH,Number=1,Type=Integer,Description=\"Reads classified as methylated\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=UNMETH,Number=1,Type=Integer,Description=\"Reads classified as unmethylated\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=FRAC,Number=1,Type=Float,Description=\"Methylated read fraction\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=MEANP,Number=1,Type=Float,Description=\"Mean probability across methylated calls\">"
    )?;
    writeln!(
        writer,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
    )?;
    Ok(())
}

fn format_alt(label: &str) -> String {
    let sanitized: String = label
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    format!("<{sanitized}>")
}

struct OutputPaths {
    write_path: PathBuf,
    compress: bool,
}

fn resolve_output_paths(path: PathBuf, compress: bool) -> OutputPaths {
    if !compress {
        return OutputPaths {
            write_path: path,
            compress: false,
        };
    }

    let is_gz = path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    if is_gz {
        OutputPaths {
            write_path: strip_gz(&path),
            compress: true,
        }
    } else {
        OutputPaths {
            write_path: path,
            compress: true,
        }
    }
}

fn strip_gz(path: &Path) -> PathBuf {
    let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
        return path.to_path_buf();
    };
    if let Some(stripped) = name.strip_suffix(".gz") {
        path.with_file_name(stripped)
    } else {
        path.to_path_buf()
    }
}

fn append_gz(path: &Path) -> PathBuf {
    let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
        return path.with_extension("gz");
    };
    let gz_name = format!("{name}.gz");
    path.with_file_name(gz_name)
}

fn bgzip_file(input: &Path) -> Result<PathBuf> {
    let status = Command::new("bgzip")
        .arg("-f")
        .arg(input)
        .status()
        .with_context(|| "failed to execute bgzip (install htslib or pass --plain-output)")?;
    if !status.success() {
        bail!("bgzip exited with {}", status);
    }

    let gz_path = append_gz(input);
    if !gz_path.exists() {
        bail!("bgzip completed but {} was not created", gz_path.display());
    }
    Ok(gz_path)
}

fn tabix_vcf(path: &Path) -> Result<()> {
    let status = Command::new("tabix")
        .args(["-f", "-p", "vcf"])
        .arg(path)
        .status()
        .with_context(|| "failed to execute tabix (install htslib or pass --plain-output)")?;
    if !status.success() {
        bail!("tabix exited with {}", status);
    }
    Ok(())
}

struct SiteTable {
    sites: BTreeMap<SiteKey, SiteStats>,
}

impl SiteTable {
    fn new() -> Self {
        Self {
            sites: BTreeMap::new(),
        }
    }

    fn record_call(
        &mut self,
        contig: &str,
        position: i64,
        motif_idx: usize,
        methylated: bool,
        probability: Option<f32>,
    ) {
        let key = SiteKey {
            contig: contig.to_string(),
            position,
            motif_idx,
        };
        let stats = self.sites.entry(key).or_insert_with(SiteStats::default);
        if methylated {
            stats.methylated += 1;
            if let Some(prob) = probability {
                stats.methylated_prob_sum += f64::from(prob);
                stats.methylated_prob_count += 1;
            }
        } else {
            stats.unmethylated += 1;
        }
    }

    fn iter(&self) -> impl Iterator<Item = (&SiteKey, &SiteStats)> {
        self.sites.iter()
    }

    fn contigs(&self) -> BTreeSet<String> {
        self.sites.keys().map(|k| k.contig.clone()).collect()
    }

    fn len(&self) -> usize {
        self.sites.len()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct SiteKey {
    contig: String,
    position: i64,
    motif_idx: usize,
}

#[derive(Default)]
struct SiteStats {
    methylated: usize,
    unmethylated: usize,
    methylated_prob_sum: f64,
    methylated_prob_count: usize,
}

impl SiteStats {
    fn depth(&self) -> usize {
        self.methylated + self.unmethylated
    }

    fn methylation_fraction(&self) -> f64 {
        let depth = self.depth();
        if depth == 0 {
            0.0
        } else {
            self.methylated as f64 / depth as f64
        }
    }

    fn mean_probability(&self) -> Option<f64> {
        (self.methylated_prob_count > 0)
            .then_some(self.methylated_prob_sum / self.methylated_prob_count as f64)
    }
}
