use std::fmt;

/// Collection of base modification calls associated with a sequencing read.
#[derive(Debug, Clone)]
pub struct ReadRecord {
    pub id: String,
    pub contig: Option<String>,
    pub strand: Option<Strand>,
    pub sequence: Vec<u8>,
    pub quality_scores: Vec<u8>,
    pub modifications: Vec<ModificationCall>,
    reference_positions: Vec<Option<i64>>,
}

impl ReadRecord {
    pub fn new(
        id: String,
        contig: Option<String>,
        strand: Option<Strand>,
        sequence: Vec<u8>,
        modifications: Vec<ModificationCall>,
        quality_scores: Vec<u8>,
        reference_positions: Vec<Option<i64>>,
    ) -> Self {
        Self {
            id,
            contig,
            strand,
            sequence,
            quality_scores,
            modifications,
            reference_positions,
        }
    }

    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    pub fn quality_scores(&self) -> &[u8] {
        &self.quality_scores
    }

    pub fn contig_name(&self) -> Option<&str> {
        self.contig.as_deref()
    }

    pub fn modification_at<'a>(
        &'a self,
        position: usize,
        label: &str,
    ) -> Option<&'a ModificationCall> {
        self.modifications
            .iter()
            .find(|call| call.position == position && call.kind.label == label)
    }

    pub fn reference_position(&self, read_position: usize) -> Option<i64> {
        self.reference_positions.get(read_position).and_then(|p| *p)
    }
}

#[derive(Debug, Clone)]
pub struct ModificationCall {
    pub position: usize,
    pub probability: Option<f32>,
    pub methylated: bool,
    pub kind: ModificationKind,
}

#[derive(Debug, Clone)]
pub struct ModificationKind {
    pub canonical_base: u8,
    pub label: String,
    pub code: String,
    pub strand: Strand,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}
