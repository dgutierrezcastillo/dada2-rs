use std::path::Path;
use anyhow::{anyhow, Result};
use needletail::parse_fastx_file;
use crate::models::DadaRecord;

/// High-performance FASTQ/FASTA sequence parser.
/// Parses biological records and returns a vector of DadaRecords.
/// Automatically handles gzip compression.
pub fn parse_sequence_file<P: AsRef<Path>>(path: P) -> Result<Vec<DadaRecord>> {
    let mut records = Vec::new();
    let mut reader = parse_fastx_file(path.as_ref())
        .map_err(|e| anyhow!("Failed to parse sequence file: {}", e))?;

    while let Some(record) = reader.next() {
        let seq_record = record.map_err(|e| anyhow!("Failed to read record: {}", e))?;
        
        // Parse ID (header)
        let id = String::from_utf8_lossy(seq_record.id()).trim().to_string();
        
        // Parse sequence and normalize to uppercase ASCII bytes
        let mut sequence = seq_record.seq().to_vec();
        sequence.make_ascii_uppercase();
        
        // Parse quality scores. If FASTA, assign dummy Phred 40 (standard high quality)
        let qualities = match seq_record.qual() {
            Some(qual) => qual.to_vec(),
            None => vec![40u8; sequence.len()],
        };

        if sequence.len() != qualities.len() {
            return Err(anyhow!(
                "Sequence length ({}) does not match quality scores length ({}) for record {}",
                sequence.len(),
                qualities.len(),
                id
            ));
        }

        records.push(DadaRecord {
            id,
            sequence,
            qualities,
        });
    }

    Ok(records)
}
