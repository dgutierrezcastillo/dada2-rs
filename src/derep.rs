use std::collections::HashMap;
use anyhow::{anyhow, Result};
use crate::models::{DadaRecord, Derep, UniqueRecord};

/// Dereplicates a slice of DadaRecords by merging identical sequences.
/// Computes abundance and consensus positional quality scores.
/// Returns a Derep object with uniques sorted by abundance in descending order.
pub fn dereplicate(records: &[DadaRecord]) -> Result<Derep> {
    if records.is_empty() {
        return Ok(Derep {
            uniques: Vec::new(),
            total_reads: 0,
        });
    }

    // Map: Sequence -> (Abundance, Accumulator for positional quality scores)
    let mut seq_map: HashMap<Vec<u8>, (u32, Vec<u64>)> = HashMap::new();

    for record in records {
        let seq = &record.sequence;
        let quals = &record.qualities;
        
        let entry = seq_map.entry(seq.clone()).or_insert_with(|| {
            (0, vec![0u64; seq.len()])
        });
        
        // Increment abundance
        entry.0 += 1;
        
        // Accumulate qualities. Ensure sequence length matches qualities length
        if entry.1.len() != quals.len() {
            return Err(anyhow!(
                "Sequence length mismatch during dereplication. Expected {}, found {} for read ID {}",
                entry.1.len(),
                quals.len(),
                record.id
            ));
        }
        
        for (acc, &q) in entry.1.iter_mut().zip(quals.iter()) {
            *acc += q as u64;
        }
    }

    // Convert the hash map into a vector of UniqueRecords
    let mut uniques = Vec::with_capacity(seq_map.len());
    for (seq, (abundance, quality_acc)) in seq_map {
        let mut average_qualities = Vec::with_capacity(quality_acc.len());
        let mut consensus_qualities = Vec::with_capacity(quality_acc.len());
        
        for acc in quality_acc {
            let avg = acc as f64 / abundance as f64;
            average_qualities.push(avg);
            consensus_qualities.push(avg.round() as u8);
        }

        uniques.push(UniqueRecord {
            sequence: seq,
            abundance,
            average_qualities,
            consensus_qualities,
        });
    }

    // Sort uniques by abundance in descending order
    uniques.sort_by(|a, b| b.abundance.cmp(&a.abundance));

    let total_reads = records.len() as u32;

    Ok(Derep {
        uniques,
        total_reads,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dereplicate_empty() {
        let result = dereplicate(&[]).unwrap();
        assert_eq!(result.total_reads, 0);
        assert!(result.uniques.is_empty());
    }

    #[test]
    fn test_dereplicate_basic() {
        let r1 = DadaRecord {
            id: "read1".to_string(),
            sequence: b"ACTG".to_vec(),
            qualities: vec![30, 40, 20, 10],
        };
        let r2 = DadaRecord {
            id: "read2".to_string(),
            sequence: b"ACTG".to_vec(),
            qualities: vec![10, 20, 30, 40],
        };
        let r3 = DadaRecord {
            id: "read3".to_string(),
            sequence: b"AAAA".to_vec(),
            qualities: vec![40, 40, 40, 40],
        };

        let result = dereplicate(&[r1, r2, r3]).unwrap();
        assert_eq!(result.total_reads, 3);
        assert_eq!(result.uniques.len(), 2);

        // ACTG should have abundance 2 (since r1 and r2 are identical)
        // AAAA should have abundance 1
        // uniques should be sorted by abundance descending
        let u1 = &result.uniques[0]; // should be ACTG
        let u2 = &result.uniques[1]; // should be AAAA

        assert_eq!(u1.sequence, b"ACTG");
        assert_eq!(u1.abundance, 2);
        // Average qualities for ACTG: [(30+10)/2, (40+20)/2, (20+30)/2, (10+40)/2] = [20.0, 30.0, 25.0, 25.0]
        assert_eq!(u1.average_qualities, vec![20.0, 30.0, 25.0, 25.0]);
        assert_eq!(u1.consensus_qualities, vec![20, 30, 25, 25]);

        assert_eq!(u2.sequence, b"AAAA");
        assert_eq!(u2.abundance, 1);
        assert_eq!(u2.average_qualities, vec![40.0, 40.0, 40.0, 40.0]);
        assert_eq!(u2.consensus_qualities, vec![40, 40, 40, 40]);
    }
}

