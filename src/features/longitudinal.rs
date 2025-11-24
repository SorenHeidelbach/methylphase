use std::collections::{BTreeMap, HashMap};

#[derive(Debug)]
pub struct LongitudinalAnalyzer {
    block_size: i64,
    contigs: BTreeMap<String, BTreeMap<i64, BlockSummary>>,
}

impl LongitudinalAnalyzer {
    pub fn new(block_size: i64) -> Self {
        assert!(block_size > 0, "block size must be positive");
        Self {
            block_size,
            contigs: BTreeMap::new(),
        }
    }

    pub fn record_call(
        &mut self,
        contig: &str,
        position: i64,
        motif: &str,
        read_id: &str,
        methylated: bool,
    ) {
        let blocks = self.contigs.entry(contig.to_string()).or_default();
        let block_idx = position.div_euclid(self.block_size);
        let block_start = block_idx * self.block_size;
        let block_end = block_start + self.block_size;
        let summary = blocks
            .entry(block_idx)
            .or_insert_with(|| BlockSummary::new(contig.to_string(), block_start, block_end));
        summary.record_call(position, motif, read_id, methylated);
    }

    pub fn finish(self) -> Vec<BlockResult> {
        let mut results = Vec::new();
        for (_, blocks) in self.contigs {
            for (_, summary) in blocks {
                results.push(summary.into_result());
            }
        }
        results
    }
}

#[derive(Debug)]
pub struct BlockResult {
    pub contig: String,
    pub block_start: i64,
    pub block_end: i64,
    pub sites: Vec<BlockSite>,
    pub patterns: Vec<PatternSummary>,
    pub statistics: BlockStatistics,
}

#[derive(Debug)]
pub struct BlockStatistics {
    pub total_reads: usize,
    pub informative_reads: usize,
    pub uninformative_reads: usize,
    pub fully_methylated_reads: usize,
    pub fully_unmethylated_reads: usize,
    pub mixed_reads: usize,
    pub mean_methylation: Option<f64>,
    pub methylation_variance: Option<f64>,
    pub polarization: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct BlockSite {
    pub position: i64,
    pub motif: String,
}

impl PartialEq for BlockSite {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position && self.motif == other.motif
    }
}

impl Eq for BlockSite {}

impl PartialOrd for BlockSite {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for BlockSite {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| self.motif.cmp(&other.motif))
    }
}

#[derive(Debug)]
pub struct PatternSummary {
    pub pattern: String,
    pub read_count: usize,
}

#[derive(Debug)]
struct BlockSummary {
    contig: String,
    block_start: i64,
    block_end: i64,
    sites: Vec<BlockSite>,
    read_patterns: HashMap<String, Vec<Option<bool>>>,
}

impl BlockSummary {
    fn new(contig: String, block_start: i64, block_end: i64) -> Self {
        Self {
            contig,
            block_start,
            block_end,
            sites: Vec::new(),
            read_patterns: HashMap::new(),
        }
    }

    fn record_call(&mut self, position: i64, motif: &str, read_id: &str, methylated: bool) {
        let site_idx = self.ensure_site(position, motif);
        let site_count = self.sites.len();
        let entry = self
            .read_patterns
            .entry(read_id.to_string())
            .or_insert_with(|| vec![None; site_count]);
        if entry.len() < site_count {
            entry.resize(site_count, None);
        }
        entry[site_idx] = Some(methylated);
    }

    fn ensure_site(&mut self, position: i64, motif: &str) -> usize {
        let key = BlockSite {
            position,
            motif: motif.to_string(),
        };
        match self.sites.binary_search(&key) {
            Ok(idx) => idx,
            Err(idx) => {
                self.sites.insert(idx, key);
                for pattern in self.read_patterns.values_mut() {
                    pattern.insert(idx, None);
                }
                idx
            }
        }
    }

    fn into_result(self) -> BlockResult {
        let total_reads = self.read_patterns.len();
        let mut informative_reads = 0usize;
        let mut fully_methylated_reads = 0usize;
        let mut fully_unmethylated_reads = 0usize;
        let mut mixed_reads = 0usize;
        let mut fractions = Vec::new();

        let mut counts: HashMap<String, usize> = HashMap::new();
        for mut pattern in self.read_patterns.into_values() {
            if pattern.len() < self.sites.len() {
                pattern.resize(self.sites.len(), None);
            }

            let mut methylated = 0usize;
            let mut unmethylated = 0usize;
            for state in &pattern {
                match state {
                    Some(true) => methylated += 1,
                    Some(false) => unmethylated += 1,
                    None => {}
                }
            }

            if methylated + unmethylated > 0 {
                informative_reads += 1;
                if unmethylated == 0 {
                    fully_methylated_reads += 1;
                } else if methylated == 0 {
                    fully_unmethylated_reads += 1;
                } else {
                    mixed_reads += 1;
                }
                let fraction = methylated as f64 / (methylated + unmethylated) as f64;
                fractions.push(fraction);
            }

            let key: String = pattern
                .iter()
                .map(|state| match state {
                    Some(true) => 'M',
                    Some(false) => 'U',
                    None => '.',
                })
                .collect();
            *counts.entry(key).or_insert(0) += 1;
        }

        let uninformative_reads = total_reads.saturating_sub(informative_reads);

        let stats = if informative_reads > 0 {
            let mean = fractions.iter().sum::<f64>() / informative_reads as f64;
            let variance = fractions
                .iter()
                .map(|f| {
                    let diff = f - mean;
                    diff * diff
                })
                .sum::<f64>()
                / informative_reads as f64;
            let polarization = fractions
                .iter()
                .map(|f| ((f - 0.5).abs() * 2.0).min(1.0))
                .sum::<f64>()
                / informative_reads as f64;

            BlockStatistics {
                total_reads,
                informative_reads,
                uninformative_reads,
                fully_methylated_reads,
                fully_unmethylated_reads,
                mixed_reads,
                mean_methylation: Some(mean),
                methylation_variance: Some(variance),
                polarization: Some(polarization),
            }
        } else {
            BlockStatistics {
                total_reads,
                informative_reads,
                uninformative_reads,
                fully_methylated_reads,
                fully_unmethylated_reads,
                mixed_reads,
                mean_methylation: None,
                methylation_variance: None,
                polarization: None,
            }
        };

        let mut patterns: Vec<PatternSummary> = counts
            .into_iter()
            .map(|(pattern, read_count)| PatternSummary {
                pattern,
                read_count,
            })
            .collect();
        patterns.sort_by(|a, b| {
            b.read_count
                .cmp(&a.read_count)
                .then_with(|| a.pattern.cmp(&b.pattern))
        });

        BlockResult {
            contig: self.contig,
            block_start: self.block_start,
            block_end: self.block_end,
            sites: self.sites,
            patterns,
            statistics: stats,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn aggregates_read_patterns_per_block() {
        let mut analyzer = LongitudinalAnalyzer::new(1000);
        analyzer.record_call("ctg", 10, "GATC_6mA_1", "read1", true);
        analyzer.record_call("ctg", 20, "GATC_6mA_1", "read1", true);
        analyzer.record_call("ctg", 10, "GATC_6mA_1", "read2", false);
        analyzer.record_call("ctg", 20, "GATC_6mA_1", "read2", false);
        analyzer.record_call("ctg", 20, "GATC_6mA_1", "read3", true);

        let results = analyzer.finish();
        assert_eq!(results.len(), 1);
        let block = &results[0];
        assert_eq!(block.sites.len(), 2);
        assert_eq!(block.patterns.len(), 3);
        let patterns: HashMap<_, _> = block
            .patterns
            .iter()
            .map(|p| (p.pattern.clone(), p.read_count))
            .collect();
        assert_eq!(patterns.get("MM"), Some(&1));
        assert_eq!(patterns.get("UU"), Some(&1));
        assert_eq!(patterns.get(".M"), Some(&1));

        assert_eq!(block.statistics.total_reads, 3);
        assert_eq!(block.statistics.informative_reads, 3);
        assert_eq!(block.statistics.fully_methylated_reads, 2);
        assert_eq!(block.statistics.fully_unmethylated_reads, 1);
        assert_eq!(block.statistics.mixed_reads, 0);
        assert!(block
            .statistics
            .mean_methylation
            .map(|m| (m - (2.0 / 3.0)).abs() < 1e-6)
            .unwrap_or(false));
        assert!(block
            .statistics
            .polarization
            .map(|p| (p - 1.0).abs() < 1e-6)
            .unwrap_or(false));
    }

    #[test]
    fn polarization_reflects_uniform_reads() {
        let mut analyzer = LongitudinalAnalyzer::new(1000);
        analyzer.record_call("ctg", 10, "GATC_6mA_1", "read1", true);
        analyzer.record_call("ctg", 20, "GATC_6mA_1", "read1", false);
        analyzer.record_call("ctg", 10, "GATC_6mA_1", "read2", false);
        analyzer.record_call("ctg", 20, "GATC_6mA_1", "read2", true);

        let results = analyzer.finish();
        let block = &results[0];
        assert_eq!(block.statistics.polarization, Some(0.0));
        assert_eq!(block.statistics.mixed_reads, 2);
    }
}
