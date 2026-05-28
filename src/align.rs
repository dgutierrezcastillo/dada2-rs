/// Represents the outcome of an ends-free sequence alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignResult {
    /// The optimal alignment score.
    pub score: i32,
    /// The aligned query sequence (with gaps represented as `-`).
    pub query_align: String,
    /// The aligned reference sequence (with gaps represented as `-`).
    pub ref_align: String,
}

/// Identifies homopolymer locations in a sequence.
/// A position is part of a homopolymer if it is within a run of at least 3 identical bases (e.g. "AAA").
/// Returns a boolean vector of the same length as the sequence.
pub fn find_homopolymers(seq: &[u8]) -> Vec<bool> {
    let len = seq.len();
    let mut homo = vec![false; len];
    if len < 3 {
        return homo;
    }

    let mut start = 0;
    for i in 0..len {
        if i == len - 1 || seq[i] != seq[i + 1] {
            let run_len = i - start + 1;
            if run_len >= 3 {
                for k in start..=i {
                    homo[k] = true;
                }
            }
            start = i + 1;
        }
    }
    
    homo
}

/// Implements a banded ends-free Needleman-Wunsch alignment.
/// Start and end gaps are not penalized.
/// Restrictions to diagonal band of size `band` are applied.
pub fn nwalign_endsfree(
    s1: &[u8],
    s2: &[u8],
    match_score: i16,
    mismatch_score: i16,
    gap_penalty: i16,
    band: isize,
) -> AlignResult {
    let len1 = s1.len();
    let len2 = s2.len();
    let nrow = len1 + 1;
    let ncol = len2 + 1;

    let mut d = vec![0i32; nrow * ncol];
    let mut p = vec![0u8; nrow * ncol];

    // Initialize left column (up moves, path = 3)
    for i in 0..=len1 {
        d[i * ncol] = 0; // ends-free gap
        p[i * ncol] = 3;
    }

    // Initialize top row (left moves, path = 2)
    for j in 0..=len2 {
        d[j] = 0; // ends-free gap
        p[j] = 2;
    }

    // Calculate left/right bands
    let lband = if len2 > len1 {
        band
    } else if len1 > len2 {
        band + (len1 - len2) as isize
    } else {
        band
    };

    let rband = if len2 > len1 {
        band + (len2 - len1) as isize
    } else if len1 > len2 {
        band
    } else {
        band
    };

    // Apply band boundary initializations to very low scores
    if band >= 0 && (band < len1 as isize || band < len2 as isize) {
        for i in 0..=len1 {
            let idx_left = i as isize - lband - 1;
            if idx_left >= 0 {
                d[i * ncol + idx_left as usize] = -9999;
            }
            let idx_right = i as isize + rband + 1;
            if idx_right <= len2 as isize {
                d[i * ncol + idx_right as usize] = -9999;
            }
        }
    }

    // DP Main Loop
    for i in 1..=len1 {
        let (l, r) = if band >= 0 {
            let mut left_idx = i as isize - lband;
            if left_idx < 1 {
                left_idx = 1;
            }
            let mut right_idx = i as isize + rband;
            if right_idx > len2 as isize {
                right_idx = len2 as isize;
            }
            (left_idx as usize, right_idx as usize)
        } else {
            (1, len2)
        };

        for j in l..=r {
            // Score for left move (gap in s1/up sequence, s2/left sequence moves)
            let left = if i == len1 {
                d[i * ncol + j - 1] // ends-free gap
            } else {
                d[i * ncol + j - 1] + gap_penalty as i32
            };

            // Score for up move (gap in s2/left sequence, s1/up sequence moves)
            let up = if j == len2 {
                d[(i - 1) * ncol + j] // ends-free gap
            } else {
                d[(i - 1) * ncol + j] + gap_penalty as i32
            };

            // Score for diagonal move
            let score_diag = if s1[i - 1] == s2[j - 1] {
                match_score as i32
            } else {
                mismatch_score as i32
            };
            let diag = d[(i - 1) * ncol + j - 1] + score_diag;

            // Tie breaking: up > left > diag
            if up >= diag && up >= left {
                d[i * ncol + j] = up;
                p[i * ncol + j] = 3;
            } else if left >= diag {
                d[i * ncol + j] = left;
                p[i * ncol + j] = 2;
            } else {
                d[i * ncol + j] = diag;
                p[i * ncol + j] = 1;
            }
        }
    }

    // Traceback
    let mut query_align = Vec::with_capacity(len1 + len2);
    let mut ref_align = Vec::with_capacity(len1 + len2);
    let mut i = len1;
    let mut j = len2;

    while i > 0 || j > 0 {
        match p[i * ncol + j] {
            1 => {
                i -= 1;
                j -= 1;
                query_align.push(s1[i]);
                ref_align.push(s2[j]);
            }
            2 => {
                j -= 1;
                query_align.push(b'-');
                ref_align.push(s2[j]);
            }
            3 => {
                i -= 1;
                query_align.push(s1[i]);
                ref_align.push(b'-');
            }
            _ => panic!("Needleman-Wunsch path out of range"),
        }
    }

    query_align.reverse();
    ref_align.reverse();

    let score = d[len1 * ncol + len2];

    AlignResult {
        score,
        query_align: String::from_utf8_lossy(&query_align).into_owned(),
        ref_align: String::from_utf8_lossy(&ref_align).into_owned(),
    }
}

/// Implements a banded ends-free Needleman-Wunsch alignment with support for homopolymer gap discounts.
pub fn nwalign_endsfree_homo(
    s1: &[u8],
    s2: &[u8],
    match_score: i16,
    mismatch_score: i16,
    gap_penalty: i16,
    homo_gap_penalty: i16,
    band: isize,
) -> AlignResult {
    let len1 = s1.len();
    let len2 = s2.len();
    let nrow = len1 + 1;
    let ncol = len2 + 1;

    let homo1 = find_homopolymers(s1);
    let homo2 = find_homopolymers(s2);

    let mut d = vec![0i32; nrow * ncol];
    let mut p = vec![0u8; nrow * ncol];

    // Initialize left column (up moves, path = 3)
    for i in 0..=len1 {
        d[i * ncol] = 0; // ends-free gap
        p[i * ncol] = 3;
    }

    // Initialize top row (left moves, path = 2)
    for j in 0..=len2 {
        d[j] = 0; // ends-free gap
        p[j] = 2;
    }

    // Calculate left/right bands
    let lband = if len2 > len1 {
        band
    } else if len1 > len2 {
        band + (len1 - len2) as isize
    } else {
        band
    };

    let rband = if len2 > len1 {
        band + (len2 - len1) as isize
    } else if len1 > len2 {
        band
    } else {
        band
    };

    // Apply band boundary initializations
    if band >= 0 && (band < len1 as isize || band < len2 as isize) {
        for i in 0..=len1 {
            let idx_left = i as isize - lband - 1;
            if idx_left >= 0 {
                d[i * ncol + idx_left as usize] = -9999;
            }
            let idx_right = i as isize + rband + 1;
            if idx_right <= len2 as isize {
                d[i * ncol + idx_right as usize] = -9999;
            }
        }
    }

    // DP Main Loop
    for i in 1..=len1 {
        let (l, r) = if band >= 0 {
            let mut left_idx = i as isize - lband;
            if left_idx < 1 {
                left_idx = 1;
            }
            let mut right_idx = i as isize + rband;
            if right_idx > len2 as isize {
                right_idx = len2 as isize;
            }
            (left_idx as usize, right_idx as usize)
        } else {
            (1, len2)
        };

        for j in l..=r {
            // Score for left move
            let left = if i == len1 {
                d[i * ncol + j - 1] // ends-free gap
            } else if homo2[j - 1] {
                d[i * ncol + j - 1] + homo_gap_penalty as i32 // homopolymer gap discount
            } else {
                d[i * ncol + j - 1] + gap_penalty as i32
            };

            // Score for up move
            let up = if j == len2 {
                d[(i - 1) * ncol + j] // ends-free gap
            } else if homo1[i - 1] {
                d[(i - 1) * ncol + j] + homo_gap_penalty as i32 // homopolymer gap discount
            } else {
                d[(i - 1) * ncol + j] + gap_penalty as i32
            };

            // Score for diagonal move
            let score_diag = if s1[i - 1] == s2[j - 1] {
                match_score as i32
            } else {
                mismatch_score as i32
            };
            let diag = d[(i - 1) * ncol + j - 1] + score_diag;

            // Tie breaking: up > left > diag
            if up >= diag && up >= left {
                d[i * ncol + j] = up;
                p[i * ncol + j] = 3;
            } else if left >= diag {
                d[i * ncol + j] = left;
                p[i * ncol + j] = 2;
            } else {
                d[i * ncol + j] = diag;
                p[i * ncol + j] = 1;
            }
        }
    }

    // Traceback
    let mut query_align = Vec::with_capacity(len1 + len2);
    let mut ref_align = Vec::with_capacity(len1 + len2);
    let mut i = len1;
    let mut j = len2;

    while i > 0 || j > 0 {
        match p[i * ncol + j] {
            1 => {
                i -= 1;
                j -= 1;
                query_align.push(s1[i]);
                ref_align.push(s2[j]);
            }
            2 => {
                j -= 1;
                query_align.push(b'-');
                ref_align.push(s2[j]);
            }
            3 => {
                i -= 1;
                query_align.push(s1[i]);
                ref_align.push(b'-');
            }
            _ => panic!("Needleman-Wunsch path out of range"),
        }
    }

    query_align.reverse();
    ref_align.reverse();

    let score = d[len1 * ncol + len2];

    AlignResult {
        score,
        query_align: String::from_utf8_lossy(&query_align).into_owned(),
        ref_align: String::from_utf8_lossy(&ref_align).into_owned(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_homopolymer_detection() {
        let seq = b"AAACCCGGT";
        let homo = find_homopolymers(seq);
        // "AAA" and "CCC" should be detected as homopolymers (positions 0-2 and 3-5)
        // "GG" should not (homopolymers must have length >= 3)
        assert_eq!(homo, vec![true, true, true, true, true, true, false, false, false]);
    }

    #[test]
    fn test_nwalign_endsfree_basic() {
        let s1 = b"ACTG";
        let s2 = b"ACTG";
        // Perfect match
        let res = nwalign_endsfree(s1, s2, 5, -4, -8, -1);
        assert_eq!(res.score, 20); // 4 * 5 = 20
        assert_eq!(res.query_align, "ACTG");
        assert_eq!(res.ref_align, "ACTG");
    }

    #[test]
    fn test_nwalign_endsfree_gaps() {
        // Ends free means start/end gaps are not penalized
        let s1 = b"ACTG";
        let s2 = b"AAACTG"; // AAA added at the beginning of s2
        let res = nwalign_endsfree(s1, s2, 5, -4, -8, -1);
        assert_eq!(res.score, 20); // Gaps on left end of s1 should be ends-free (0 penalty)
        assert_eq!(res.query_align, "--ACTG");
        assert_eq!(res.ref_align, "AAACTG");
    }

    #[test]
    fn test_nwalign_endsfree_homo() {
        let s1 = b"CTGAAAAACTG";
        let s2 = b"CTGAAACTG"; // AA deletion within homopolymer of A's in s2 (middle of sequence)
        // If we use standard alignment, gap penalty is -8
        let res_std = nwalign_endsfree(s1, s2, 5, -4, -8, -1);
        
        // If we use homopolymer alignment, homopolymer gap penalty is -2 (much lower discount)
        let res_homo = nwalign_endsfree_homo(s1, s2, 5, -4, -8, -2, -1);
        
        // Homopolymer score must be significantly better than standard score due to discount
        assert!(res_homo.score > res_std.score);
    }

}
