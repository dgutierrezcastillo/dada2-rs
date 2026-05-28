/// Represents an individual sequence read parsed from a FASTQ or FASTA file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DadaRecord {
    /// The read identifier or header.
    pub id: String,
    /// The biological sequence represented as ASCII bytes (e.g., b"ACTG").
    pub sequence: Vec<u8>,
    /// The positional numerical Phred quality scores.
    pub qualities: Vec<u8>,
}

/// Represents a unique biological sequence resolved during dereplication.
#[derive(Debug, Clone, PartialEq)]
pub struct UniqueRecord {
    /// The unique sequence represented as ASCII bytes (e.g., b"ACTG").
    pub sequence: Vec<u8>,
    /// The number of reads corresponding to this unique sequence.
    pub abundance: u32,
    /// The average floating-point quality score at each sequence position.
    pub average_qualities: Vec<f64>,
    /// The rounded integer Phred quality score at each sequence position (matching C++ Raw.qual).
    pub consensus_qualities: Vec<u8>,
}

/// Container representing a dereplicated sample, matching R's `derep-class` output.
#[derive(Debug, Clone, PartialEq)]
pub struct Derep {
    /// List of all unique sequence records sorted by abundance in descending order.
    pub uniques: Vec<UniqueRecord>,
    /// Total number of reads parsed from the original file.
    pub total_reads: u32,
}
