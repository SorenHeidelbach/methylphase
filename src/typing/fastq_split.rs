use crate::typing::error::MethylError;
use regex::Regex;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Split a FASTQ into per-class FASTQ files based on assignments TSV.
///
/// Assignments file must have an `id` column and either:
/// - `assigned_class` column, or
/// - responsibility columns named `class_0`, `resp_class_0`, etc. (argmax used).
pub fn split_fastq(
    fastq: &Path,
    assignments: &Path,
    delimiter: u8,
    out_dir: &Path,
) -> Result<(), MethylError> {
    fs::create_dir_all(out_dir)?;
    let assign_map = load_assignments(assignments, delimiter)?;

    let fh = File::open(fastq)?;
    let reader = BufReader::new(fh);

    let mut writers: HashMap<usize, BufWriter<File>> = HashMap::new();
    let mut unknown = BufWriter::new(File::create(out_dir.join("class_unknown.fastq"))?);

    let mut lines = reader.lines();
    while let Some(h) = lines.next() {
        let header = h?;
        let seq = lines.next().ok_or_else(|| MethylError::Parse("truncated FASTQ".into()))??;
        let plus = lines.next().ok_or_else(|| MethylError::Parse("truncated FASTQ".into()))??;
        let qual = lines.next().ok_or_else(|| MethylError::Parse("truncated FASTQ".into()))??;

        let rid = parse_read_id(&header);
        if let Some(cls) = assign_map.get(&rid) {
            let writer = writers.entry(*cls).or_insert_with(|| {
                let path = out_dir.join(format!("class_{}.fastq", cls));
                BufWriter::new(File::create(path).expect("create fastq"))
            });
            write_record(writer, &header, &seq, &plus, &qual)?;
        } else {
            write_record(&mut unknown, &header, &seq, &plus, &qual)?;
        }
    }
    unknown.flush()?;
    for (_, mut w) in writers.into_iter() {
        w.flush()?;
    }
    Ok(())
}

fn write_record<W: Write>(
    w: &mut W,
    h: &str,
    s: &str,
    p: &str,
    q: &str,
) -> Result<(), MethylError> {
    writeln!(w, "{}", h)?;
    writeln!(w, "{}", s)?;
    writeln!(w, "{}", p)?;
    writeln!(w, "{}", q)?;
    Ok(())
}

fn parse_read_id(header: &str) -> String {
    // Use token after '@' up to first whitespace.
    header
        .trim_start_matches('@')
        .split_whitespace()
        .next()
        .unwrap_or(header)
        .to_string()
}

fn load_assignments(
    path: &Path,
    delimiter: u8,
) -> Result<HashMap<String, usize>, MethylError> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(true)
        .from_path(path)?;

    let headers = reader.headers()?.clone();
    let id_idx = headers
        .iter()
        .position(|h| h == "id")
        .ok_or_else(|| MethylError::Parse("assignments missing id column".into()))?;
    let assigned_idx = headers.iter().position(|h| h == "assigned_class");

    // Responsibility columns
    let resp_re = Regex::new(r"^(class_|resp_class_)(\d+)$").unwrap();
    let mut resp_cols: Vec<(usize, usize)> = Vec::new(); // (col_idx, class)
    if assigned_idx.is_none() {
        for (i, h) in headers.iter().enumerate() {
            if let Some(caps) = resp_re.captures(h) {
                if let Some(m) = caps.get(2) {
                    if let Ok(cls) = m.as_str().parse::<usize>() {
                        resp_cols.push((i, cls));
                    }
                }
            }
        }
    }

    if assigned_idx.is_none() && resp_cols.is_empty() {
        return Err(MethylError::Parse(
            "assignments file must have assigned_class or class_*/resp_class_* columns".into(),
        ));
    }

    let mut map = HashMap::new();
    for rec in reader.records() {
        let rec = rec?;
        let id = rec
            .get(id_idx)
            .ok_or_else(|| MethylError::Parse("missing id value".into()))?
            .to_string();
        let cls = if let Some(idx) = assigned_idx {
            let c: usize = rec
                .get(idx)
                .ok_or_else(|| MethylError::Parse("missing assigned_class".into()))?
                .parse()
                .map_err(|e| MethylError::Parse(format!("invalid assigned_class: {}", e)))?;
            c
        } else {
            // argmax responsibility
            let mut best_cls = 0;
            let mut best_val = -1.0;
            for (col, cls) in &resp_cols {
                if let Some(vs) = rec.get(*col) {
                    if let Ok(v) = vs.parse::<f64>() {
                        if v > best_val {
                            best_val = v;
                            best_cls = *cls;
                        }
                    }
                }
            }
            best_cls
        };
        map.insert(id, cls);
    }
    Ok(map)
}
