use anyhow::{Context, Result};
use bstr::ByteSlice;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_sam::{self as sam};
use std::{borrow::Cow, fs::File, path::Path};

/// Thin wrapper around an indexed BAM reader that exposes raw records.
pub struct BamReader {
    reader: bam::io::IndexedReader<bgzf::io::Reader<File>>,
    header: sam::Header,
    reference_names: Vec<String>,
}

impl BamReader {
    pub fn open<P>(bam_path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let bam_path = bam_path.as_ref();
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .with_context(|| format!("could not open BAM {}", bam_path.display()))?;
        let header = reader.read_header().context("failed to read BAM header")?;
        let reference_names = header
            .reference_sequences()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect();

        Ok(Self {
            reader,
            header,
            reference_names,
        })
    }

    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    pub fn contigs(&self) -> Vec<String> {
        self.reference_names.clone()
    }

    pub fn read_ids(&mut self, contig: &str, limit: Option<usize>) -> Result<Vec<String>> {
        let region: Region = contig
            .parse()
            .with_context(|| format!("invalid contig or region: {contig}"))?;
        let mut query = self
            .reader
            .query(&self.header, &region)
            .with_context(|| format!("failed to query region {contig}"))?;

        let mut reads = Vec::new();
        while let Some(result) = query.next() {
            let record = result.context("BAM record read failed")?;
            if let Some(name) = record.name() {
                let raw_name: &[u8] = name.as_ref();
                let read_id = raw_name
                    .to_str()
                    .context("read name contains invalid UTF-8")?;
                reads.push(read_id.to_string());
            }

            if let Some(max) = limit {
                if reads.len() >= max {
                    break;
                }
            }
        }

        Ok(reads)
    }

    pub fn visit_records<F>(&mut self, contigs: Option<&[String]>, mut visitor: F) -> Result<()>
    where
        F: FnMut(bam::Record) -> Result<()>,
    {
        let targets: Cow<'_, [String]> = match contigs {
            Some(cs) if !cs.is_empty() => Cow::Borrowed(cs),
            _ => Cow::Owned(self.contigs()),
        };

        for contig in targets.iter() {
            self.visit_contig(contig, &mut visitor)?;
        }

        Ok(())
    }

    fn visit_contig<F>(&mut self, contig: &str, visitor: &mut F) -> Result<()>
    where
        F: FnMut(bam::Record) -> Result<()>,
    {
        let region: Region = contig
            .parse()
            .with_context(|| format!("invalid contig or region: {contig}"))?;
        let mut query = self
            .reader
            .query(&self.header, &region)
            .with_context(|| format!("failed to query region {contig}"))?;

        while let Some(result) = query.next() {
            let record = result.context("BAM record read failed")?;
            visitor(record)?;
        }

        Ok(())
    }
}
