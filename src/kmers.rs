/// Helper function to convert a nucleotide byte to a 2-bit integer.
/// Returns Some(index) for standard nucleotides, or None for invalid characters (like 'N' or gaps).
#[inline]
pub fn nucleotide_to_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Generates a k-mer count vector for a given sequence and k-mer size.
/// For standard k = 5, the vector has a size of 4^5 = 1024.
pub fn assign_kmer(seq: &[u8], k: usize) -> Vec<u16> {
    let n_kmer = 1 << (2 * k); // 4^k possible kmers
    let mut kvec = vec![0u16; n_kmer];
    if seq.len() < k {
        return kvec;
    }

    let klen = seq.len() - k + 1;
    for i in 0..klen {
        let mut kmer_idx = 0;
        let mut valid = true;
        
        for j in 0..k {
            match nucleotide_to_index(seq[i + j]) {
                Some(idx) => {
                    kmer_idx = 4 * kmer_idx + idx;
                }
                None => {
                    valid = false;
                    break;
                }
            }
        }
        
        if valid {
            kvec[kmer_idx] += 1;
        }
    }
    
    kvec
}

/// Calculates the k-mer distance between two unique sequences.
/// The distance metric matches the C++ kmer_dist implementation exactly.
pub fn kmer_dist(kv1: &[u16], len1: usize, kv2: &[u16], len2: usize, k: usize) -> f64 {
    let n_kmer = 1 << (2 * k);
    let min_len = if len1 < len2 { len1 } else { len2 };
    
    // Safety check: if min_len is less than k, the distance is maximum (1.0)
    if min_len < k {
        return 1.0;
    }

    // Accumulate the sum of the minimums of the two k-mer count vectors.
    // Highly vectorizable loop.
    let mut dotsum = 0u64;
    for i in 0..n_kmer {
        let v1 = kv1[i];
        let v2 = kv2[i];
        dotsum += if v1 < v2 { v1 } else { v2 } as u64;
    }
    
    let denominator = (min_len - k + 1) as f64;
    let dot = dotsum as f64 / denominator;
    
    1.0 - dot
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_basic() {
        let seq1 = b"AAAAA"; // Only one 5-mer (AAAAA -> index 0)
        let kv1 = assign_kmer(seq1, 5);
        assert_eq!(kv1[0], 1);
        assert_eq!(kv1.iter().sum::<u16>(), 1);

        let seq2 = b"AAAAC"; // 5-mer AAAAC -> index 1
        let kv2 = assign_kmer(seq2, 5);
        assert_eq!(kv2[1], 1);
        assert_eq!(kv2.iter().sum::<u16>(), 1);

        // Identical sequence distance must be 0.0
        let dist_self = kmer_dist(&kv1, seq1.len(), &kv1, seq1.len(), 5);
        assert_eq!(dist_self, 0.0);

        // Completely divergent distance must be 1.0
        let dist_div = kmer_dist(&kv1, seq1.len(), &kv2, seq2.len(), 5);
        assert_eq!(dist_div, 1.0);
    }

    #[test]
    fn test_kmer_overlapping() {
        let seq1 = b"AAAAAA"; // Two overlapping AAAAA (index 0)
        let kv1 = assign_kmer(seq1, 5);
        assert_eq!(kv1[0], 2);

        let seq2 = b"AAAAAC"; // One AAAAA, one AAAAC
        let kv2 = assign_kmer(seq2, 5);
        assert_eq!(kv2[0], 1);
        assert_eq!(kv2[1], 1);

        // Common count sum = min(2, 1) + min(0, 1) = 1 + 0 = 1
        // Denominator = 6 - 5 + 1 = 2
        // Dist = 1.0 - (1.0 / 2.0) = 0.5
        let dist = kmer_dist(&kv1, seq1.len(), &kv2, seq2.len(), 5);
        assert_eq!(dist, 0.5);
    }
}
