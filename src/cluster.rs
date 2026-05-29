use crate::align;
use crate::models::UniqueRecord;
use rayon::prelude::*;
use statrs::distribution::{DiscreteCDF, Poisson};

/// Represents a comparison between a cluster center and a unique sequence.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Comparison {
    /// The cluster index (C++ Bi index `i`).
    pub cluster_index: usize,
    /// The unique record index in the clustering records list.
    pub raw_index: usize,
    /// The calculated error transition rate (product of base transitions).
    pub lambda: f64,
    /// The Hamming distance (number of substitutions, skipping gaps).
    pub hamming: u32,
}

/// Represents a tracking record for a unique sequence during divisive partitioning.
#[derive(Debug, Clone)]
pub struct ClusteringRecord {
    /// The original unique record index in `Derep.uniques`.
    pub index: usize,
    /// The sequence abundance (read count).
    pub abundance: u32,
    /// The abundance p-value relative to its assigned cluster.
    pub p_value: f64,
    /// Locked to its current cluster in greedy mode.
    pub lock: bool,
    /// The comparison to its assigned cluster center.
    pub comp: Option<Comparison>,
    /// The maximum expected reads (E_minmax) across all clusters.
    pub e_minmax: f64,
    /// Prior reason to expect this sequence (e.g. from an existing reference).
    pub prior: bool,
}

/// Represents a partition or cluster (C++ `Bi` struct).
#[derive(Debug, Clone)]
pub struct Cluster {
    /// The cluster index in the clustering list.
    pub index: usize,
    /// The index of the unique sequence serving as the partition center.
    pub center_index: usize,
    /// List of unique sequence indices assigned to this cluster.
    pub raw_indices: Vec<usize>,
    /// Total reads (sum of abundances of all sequences in the cluster).
    pub total_reads: u32,
    /// Self-production error rate (lambda of center aligned to itself).
    pub self_lambda: f64,
    /// Pairwise comparisons between all unique sequences and this cluster center.
    pub comparisons: Vec<Option<Comparison>>,
    /// Flagged to recalculate expected abundances and update.
    pub update_e: bool,
    /// Flagged to check and lock raws in greedy mode.
    pub check_locks: bool,
    
    // Birth metadata for lineage tracking:
    pub birth_type: String, // "I" for initial, "A" for abundance, "P" for prior
    pub birth_from: Option<usize>,
    pub birth_pval: f64,
    pub birth_fold: f64,
    pub birth_e: f64,
}

/// Container for the entire divisive partitioning system (C++ `B` struct).
#[derive(Debug, Clone)]
pub struct Clustering {
    /// Tracking records for each unique sequence.
    pub records: Vec<ClusteringRecord>,
    /// Active partitions (clusters).
    pub clusters: Vec<Cluster>,
    /// Total read count in the sample.
    pub total_reads: u32,
    /// Abundance p-value threshold (OmegaA, default 1e-40).
    pub omega_a: f64,
    /// Prior p-value threshold (OmegaP, default 1e-40).
    pub omega_p: f64,
    /// Boolean flag to use Phred quality scores.
    pub use_quals: bool,
}

/// Calculates the transition error rate lambda based on alignment and Phred quality scores.
pub fn compute_lambda(
    ref_align: &str,
    query_align: &str,
    query_qual: &[u8],
    err_mat: &ndarray::Array2<f64>,
    use_quals: bool,
) -> f64 {
    let mut lambda = 1.0;
    let mut pos1 = 0;

    let ref_bytes = ref_align.as_bytes();
    let query_bytes = query_align.as_bytes();
    let ncol = err_mat.ncols();

    for c in 0..ref_bytes.len() {
        if query_bytes[c] != b'-' {
            let q_nt = match query_bytes[c] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 0,
            };

            let r_nt = if ref_bytes[c] == b'-' {
                q_nt // Insertion: treat as match
            } else {
                match ref_bytes[c] {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => 0,
                }
            };

            let trans_idx = r_nt * 4 + q_nt;
            let q_idx = if use_quals {
                query_qual[pos1] as usize
            } else {
                0
            };

            // Clamp quality index to avoid index out of bounds
            let q_idx = q_idx.min(ncol - 1);
            lambda *= err_mat[[trans_idx, q_idx]];
            pos1 += 1;
        }
    }

    lambda.clamp(0.0, 1.0)
}

/// Generates a standard default $16 \times 41$ transition probability matrix
/// based on nominal Phred quality scores.
pub fn generate_default_error_matrix(max_q: usize) -> ndarray::Array2<f64> {
    let ncol = max_q + 1;
    let mut err_mat = ndarray::Array2::zeros((16, ncol));
    for q in 0..ncol {
        let p_err = 10.0f64.powf(-(q as f64) / 10.0);
        let p_err = p_err.min(0.75); // Cap maximum error rate
        let p_match = 1.0 - p_err;
        let p_mismatch = p_err / 3.0;

        for r in 0..4 {
            for c in 0..4 {
                let trans_idx = r * 4 + c;
                if r == c {
                    err_mat[[trans_idx, q]] = p_match;
                } else {
                    err_mat[[trans_idx, q]] = p_mismatch;
                }
            }
        }
    }
    err_mat
}

impl Clustering {
    /// Creates and initializes a new divisive partitioning session.
    pub fn new(uniques: &[UniqueRecord], omega_a: f64, omega_p: f64, use_quals: bool) -> Self {
        let nraw = uniques.len();
        let total_reads = uniques.iter().map(|u| u.abundance).sum();

        let mut records = Vec::with_capacity(nraw);
        for index in 0..nraw {
            records.push(ClusteringRecord {
                index,
                abundance: uniques[index].abundance,
                p_value: 1.0,
                lock: false,
                comp: None,
                e_minmax: 0.0,
                prior: false,
            });
        }

        let raw_indices: Vec<usize> = (0..nraw).collect();
        let initial_cluster = Cluster {
            index: 0,
            center_index: 0, // Rank 1 sequence is the most abundant center
            raw_indices,
            total_reads,
            self_lambda: 1.0,
            comparisons: vec![None; nraw],
            update_e: true,
            check_locks: true,
            birth_type: "I".to_string(),
            birth_from: None,
            birth_pval: 1.0,
            birth_fold: 1.0,
            birth_e: 0.0,
        };

        Clustering {
            records,
            clusters: vec![initial_cluster],
            total_reads,
            omega_a,
            omega_p,
            use_quals,
        }
    }

    /// Performs alignments and calculates transition rates (lambdas) for all raws to the specified cluster.
    /// This step is fully parallelized using `Rayon`.
    pub fn compare(
        &mut self,
        cluster_idx: usize,
        uniques: &[UniqueRecord],
        err_mat: &ndarray::Array2<f64>,
        match_score: i16,
        mismatch_score: i16,
        gap_penalty: i16,
        homo_gap_penalty: i16,
        band_size: isize,
        greedy: bool,
    ) {
        let center_rec_idx = self.clusters[cluster_idx].center_index;
        let center_seq = &uniques[center_rec_idx].sequence;
        let center_reads = self.records[center_rec_idx].abundance;

        // Perform parallel pairwise comparisons using Rayon
        let comps: Vec<Option<Comparison>> = (0..self.records.len())
            .into_par_iter()
            .map(|raw_idx| {
                let raw = &self.records[raw_idx];
                if greedy && (raw.abundance > center_reads) {
                    return None;
                }
                if greedy && raw.lock {
                    return None;
                }

                let raw_seq = &uniques[raw_idx].sequence;
                let align_res = align::nwalign_endsfree_homo(
                    center_seq,
                    raw_seq,
                    match_score,
                    mismatch_score,
                    gap_penalty,
                    homo_gap_penalty,
                    band_size,
                );

                let lambda = compute_lambda(
                    &align_res.query_align,
                    &align_res.ref_align,
                    &uniques[raw_idx].consensus_qualities,
                    err_mat,
                    self.use_quals,
                );

                let hamming = align_res.query_align.bytes().zip(align_res.ref_align.bytes())
                    .filter(|(q, r)| q != r && *q != b'-' && *r != b'-')
                    .count() as u32;

                Some(Comparison {
                    cluster_index: cluster_idx,
                    raw_index: raw_idx,
                    lambda,
                    hamming,
                })
            })
            .collect();

        // Update clustering records and clusters sequentially based on comparisons
        for raw_idx in 0..self.records.len() {
            if let Some(comp) = comps[raw_idx] {
                if raw_idx == center_rec_idx {
                    self.clusters[cluster_idx].self_lambda = comp.lambda;
                }

                let expected_reads_total = comp.lambda * self.total_reads as f64;
                if expected_reads_total > self.records[raw_idx].e_minmax {
                    let expected_reads_cluster = comp.lambda * center_reads as f64;
                    if expected_reads_cluster > self.records[raw_idx].e_minmax {
                        self.records[raw_idx].e_minmax = expected_reads_cluster;
                    }

                    self.clusters[cluster_idx].comparisons[raw_idx] = Some(comp);

                    if cluster_idx == 0 || raw_idx == center_rec_idx {
                        self.records[raw_idx].comp = Some(comp);
                    }
                }
            }
        }
    }

    /// Recalculates expected read rates and updates locks/p-values.
    pub fn p_update(&mut self, greedy: bool, detect_singletons: bool) {
        for i in 0..self.clusters.len() {
            let update_e = self.clusters[i].update_e;
            let center_rec_idx = self.clusters[i].center_index;
            let center_reads = self.records[center_rec_idx].abundance;

            if update_e {
                let raw_indices = self.clusters[i].raw_indices.clone();
                for &raw_idx in &raw_indices {
                    self.records[raw_idx].p_value = self.get_pa(raw_idx, i, detect_singletons);
                }
                self.clusters[i].update_e = false;
            }

            if greedy && self.clusters[i].check_locks {
                let raw_indices = self.clusters[i].raw_indices.clone();
                for &raw_idx in &raw_indices {
                    if let Some(comp) = self.clusters[i].comparisons[raw_idx] {
                        let e_reads_center = center_reads as f64 * comp.lambda;
                        if e_reads_center > self.records[raw_idx].abundance as f64 {
                            self.records[raw_idx].lock = true;
                        }
                    }
                    if raw_idx == center_rec_idx {
                        self.records[raw_idx].lock = true;
                    }
                }
                self.clusters[i].check_locks = false;
            }
        }
    }

    /// Calculates the abundance p-value for a unique sequence under the Poisson model.
    pub fn get_pa(&self, raw_idx: usize, cluster_idx: usize, detect_singletons: bool) -> f64 {
        let raw = &self.records[raw_idx];
        let bi = &self.clusters[cluster_idx];

        let comp = match bi.comparisons[raw_idx] {
            Some(c) => c,
            None => return 1.0,
        };

        let lambda = comp.lambda;
        let hamming = comp.hamming;

        if raw.abundance == 1 && !raw.prior && !detect_singletons {
            1.0
        } else if hamming == 0 {
            1.0
        } else if lambda == 0.0 {
            0.0
        } else {
            let e_reads = lambda * bi.total_reads as f64;
            self.calc_pa(raw.abundance as i32, e_reads, raw.prior || detect_singletons)
        }
    }

    /// Exact Poisson survival function normalized for conditional observation.
    pub fn calc_pa(&self, reads: i32, e_reads: f64, prior: bool) -> f64 {
        if reads <= 0 {
            return 1.0;
        }
        let n = (reads - 1) as u64;
        let dist = match Poisson::new(e_reads) {
            Ok(d) => d,
            Err(_) => return 1.0,
        };

        let cdf_val = dist.cdf(n);
        let mut pval = 1.0 - cdf_val;

        if !prior {
            let mut norm = 1.0 - (-e_reads).exp();
            if norm < 1e-7 {
                norm = e_reads - 0.5 * e_reads * e_reads;
            }
            pval = pval / norm;
        }

        pval.clamp(0.0, 1.0)
    }

    /// Reassigns each sequence to the cluster that maximizes the expected abundance.
    pub fn shuffle(&mut self) -> bool {
        let nraw = self.records.len();
        let nclust = self.clusters.len();

        let mut best_comps = vec![None; nraw];
        let mut emax = vec![0.0; nraw];

        // Initialize with cluster 0
        for index in 0..nraw {
            if let Some(comp) = self.clusters[0].comparisons[index] {
                best_comps[index] = Some(comp);
                emax[index] = comp.lambda * self.clusters[0].total_reads as f64;
            }
        }

        // Iterate over other clusters
        for i in 1..nclust {
            let bi = &self.clusters[i];
            for index in 0..nraw {
                if let Some(comp) = bi.comparisons[index] {
                    let e = comp.lambda * bi.total_reads as f64;
                    if e > emax[index] {
                        best_comps[index] = Some(comp);
                        emax[index] = e;
                    }
                }
            }
        }

        let mut shuffled = false;

        // Shuffle sequences between clusters
        for i in 0..nclust {
            let mut r = self.clusters[i].raw_indices.len() as isize - 1;
            while r >= 0 {
                let raw_idx = self.clusters[i].raw_indices[r as usize];

                if let Some(best_comp) = best_comps[raw_idx] {
                    if best_comp.cluster_index != i {
                        if raw_idx == self.clusters[i].center_index {
                            r -= 1;
                            continue;
                        }

                        // Move raw record
                        self.clusters[i].raw_indices.remove(r as usize);
                        self.clusters[best_comp.cluster_index].raw_indices.push(raw_idx);

                        self.records[raw_idx].comp = Some(best_comp);
                        shuffled = true;
                    }
                }
                r -= 1;
            }
        }

        shuffled
    }

    /// Scans p-values, finds the minimum, and divides a cluster if significant.
    pub fn bud(
        &mut self,
        min_fold: f64,
        min_hamming: u32,
        min_abund: u32,
        verbose: bool,
    ) -> Option<usize> {
        let nclust = self.clusters.len();
        let nraw = self.records.len();

        let mut min_rec_idx: Option<usize> = None;
        let mut min_cluster_idx: Option<usize> = None;
        let mut min_r_idx: Option<usize> = None;

        let mut min_rec_prior_idx: Option<usize> = None;
        let mut min_cluster_prior_idx: Option<usize> = None;
        let mut min_r_prior_idx: Option<usize> = None;

        for i in 0..nclust {
            let bi = &self.clusters[i];
            for r in 1..bi.raw_indices.len() {
                let raw_idx = bi.raw_indices[r];
                let raw = &self.records[raw_idx];

                if raw.abundance < min_abund {
                    continue;
                }

                let comp = match raw.comp {
                    Some(c) => c,
                    None => continue,
                };

                if comp.hamming >= min_hamming {
                    let expected_reads = comp.lambda * bi.total_reads as f64;
                    if min_fold <= 1.0 || raw.abundance as f64 >= min_fold * expected_reads {
                        let is_better = match min_rec_idx {
                            None => true,
                            Some(current_min_idx) => {
                                let cur_p = self.records[current_min_idx].p_value;
                                if raw.p_value < cur_p {
                                    true
                                } else if raw.p_value == cur_p && raw.abundance > self.records[current_min_idx].abundance {
                                    true
                                } else {
                                    false
                                }
                            }
                        };

                        if is_better {
                            min_rec_idx = Some(raw_idx);
                            min_cluster_idx = Some(i);
                            min_r_idx = Some(r);
                        }

                        if raw.prior {
                            let is_better_prior = match min_rec_prior_idx {
                                None => true,
                                Some(current_min_prior_idx) => {
                                    let cur_p = self.records[current_min_prior_idx].p_value;
                                    if raw.p_value < cur_p {
                                        true
                                    } else if raw.p_value == cur_p && raw.abundance > self.records[current_min_prior_idx].abundance {
                                        true
                                    } else {
                                        false
                                    }
                                }
                            };

                            if is_better_prior {
                                min_rec_prior_idx = Some(raw_idx);
                                min_cluster_prior_idx = Some(i);
                                min_r_prior_idx = Some(r);
                            }
                        }
                    }
                }
            }
        }

        // Apply Bonferroni correction and test for budding
        if let Some(raw_idx) = min_rec_idx {
            let p_a = self.records[raw_idx].p_value * nraw as f64;
            if p_a < self.omega_a {
                let mini = min_cluster_idx.unwrap();
                let minr = min_r_idx.unwrap();

                let raw_idx_popped = self.clusters[mini].raw_indices.remove(minr);
                let new_cluster_idx = self.clusters.len();
                let comp = self.records[raw_idx_popped].comp.unwrap();
                let expected = comp.lambda * self.clusters[mini].total_reads as f64;

                let new_bi = Cluster {
                    index: new_cluster_idx,
                    center_index: raw_idx_popped,
                    raw_indices: vec![raw_idx_popped],
                    total_reads: self.records[raw_idx_popped].abundance,
                    self_lambda: 1.0,
                    comparisons: vec![None; nraw],
                    update_e: true,
                    check_locks: true,
                    birth_type: "A".to_string(),
                    birth_from: Some(mini),
                    birth_pval: p_a,
                    birth_fold: self.records[raw_idx_popped].abundance as f64 / expected,
                    birth_e: expected,
                };

                self.clusters.push(new_bi);

                if verbose {
                    println!(
                        ", Division (naive): Raw {} from Bi {}, pA={:.2e}",
                        raw_idx_popped, mini, p_a
                    );
                }

                return Some(new_cluster_idx);
            }
        }

        if let Some(raw_prior_idx) = min_rec_prior_idx {
            let p_p = self.records[raw_prior_idx].p_value;
            if p_p < self.omega_p {
                let mini_prior = min_cluster_prior_idx.unwrap();
                let minr_prior = min_r_prior_idx.unwrap();

                let raw_idx_popped = self.clusters[mini_prior].raw_indices.remove(minr_prior);
                let new_cluster_idx = self.clusters.len();
                let comp = self.records[raw_idx_popped].comp.unwrap();
                let expected = comp.lambda * self.clusters[mini_prior].total_reads as f64;

                let new_bi = Cluster {
                    index: new_cluster_idx,
                    center_index: raw_idx_popped,
                    raw_indices: vec![raw_idx_popped],
                    total_reads: self.records[raw_idx_popped].abundance,
                    self_lambda: 1.0,
                    comparisons: vec![None; nraw],
                    update_e: true,
                    check_locks: true,
                    birth_type: "P".to_string(),
                    birth_from: Some(mini_prior),
                    birth_pval: p_p,
                    birth_fold: self.records[raw_idx_popped].abundance as f64 / expected,
                    birth_e: expected,
                };

                self.clusters.push(new_bi);

                if verbose {
                    println!(
                        ", Division (prior): Raw {} from Bi {}, pP={:.2e}",
                        raw_idx_popped, mini_prior, p_p
                    );
                }

                return Some(new_cluster_idx);
            }
        }

        None
    }

    /// Performs unique records census in a partition.
    pub fn bi_census(&mut self, cluster_idx: usize) {
        let mut reads = 0;
        let raw_indices = self.clusters[cluster_idx].raw_indices.clone();
        for &raw_idx in &raw_indices {
            reads += self.records[raw_idx].abundance;
        }
        if reads != self.clusters[cluster_idx].total_reads {
            self.clusters[cluster_idx].update_e = true;
        }
        self.clusters[cluster_idx].total_reads = reads;
    }

    /// Re-evaluates and assigns the cluster representative sequence center.
    pub fn bi_assign_center(&mut self, cluster_idx: usize) {
        let raw_indices = self.clusters[cluster_idx].raw_indices.clone();

        for &raw_idx in &raw_indices {
            self.records[raw_idx].lock = false;
        }

        let mut max_reads = 0;
        let mut center_rec_idx = raw_indices[0];
        for &raw_idx in &raw_indices {
            let reads = self.records[raw_idx].abundance;
            if reads > max_reads {
                max_reads = reads;
                center_rec_idx = raw_idx;
            }
        }

        self.clusters[cluster_idx].center_index = center_rec_idx;
        self.clusters[cluster_idx].check_locks = true;
    }

    /// Executes the full DADA divisive clustering pipeline.
    pub fn run_dada(
        uniques: &[UniqueRecord],
        err_mat: &ndarray::Array2<f64>,
        omega_a: f64,
        omega_p: f64,
        use_quals: bool,
        match_score: i16,
        mismatch_score: i16,
        gap_penalty: i16,
        homo_gap_penalty: i16,
        band_size: isize,
        greedy: bool,
        verbose: bool,
    ) -> Self {
        let mut clustering = Clustering::new(uniques, omega_a, omega_p, use_quals);

        if verbose {
            println!("Initializing DADA clustering with cluster 0...");
        }

        // Compare all raws to the initial center in cluster 0
        clustering.compare(
            0,
            uniques,
            err_mat,
            match_score,
            mismatch_score,
            gap_penalty,
            homo_gap_penalty,
            band_size,
            greedy,
        );

        // Update expected parameters and initial p-values
        clustering.p_update(greedy, false);

        let mut step = 0;
        loop {
            step += 1;
            if verbose {
                print!("Step {}: checking for budding...", step);
            }

            // Check if any partition center should bud
            let new_cluster_idx = clustering.bud(1.0, 1, 1, verbose);

            match new_cluster_idx {
                Some(new_i) => {
                    // Budded a new cluster: compare all records to the new center
                    clustering.compare(
                        new_i,
                        uniques,
                        err_mat,
                        match_score,
                        mismatch_score,
                        gap_penalty,
                        homo_gap_penalty,
                        band_size,
                        greedy,
                    );

                    // Reassign and shuffle raws until partitioning converges
                    loop {
                        let shuffled = clustering.shuffle();
                        if !shuffled {
                            break;
                        }

                        // Re-census and assign new centers to modified partitions
                        for i in 0..clustering.clusters.len() {
                            clustering.bi_census(i);
                            clustering.bi_assign_center(i);
                        }

                        // Re-compare partitions whose centers changed
                        for i in 0..clustering.clusters.len() {
                            if clustering.clusters[i].update_e {
                                clustering.compare(
                                    i,
                                    uniques,
                                    err_mat,
                                    match_score,
                                    mismatch_score,
                                    gap_penalty,
                                    homo_gap_penalty,
                                    band_size,
                                    greedy,
                                );
                            }
                        }
                    }

                    // Update expected rates and p-values
                    clustering.p_update(greedy, false);
                }
                None => {
                    if verbose {
                        println!("\nClustering converged with {} clusters.", clustering.clusters.len());
                    }
                    break;
                }
            }
        }

        clustering
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clustering_basic() {
        // Create 3 unique records:
        // Rank 1: 50 bases long center
        // Rank 2: 50 bases long with a substitution mismatch in the middle (A's instead of G's)
        // Rank 3: 50 bases long with a heavily mutated middle region (T's instead of G's)
        let uniques = vec![
            UniqueRecord {
                sequence: b"AAAAACCCCCGGGGGAAAAACCCCCGGGGGAAAAACCCCCGGGGGAAAAA".to_vec(),
                abundance: 500,
                average_qualities: vec![40.0; 50],
                consensus_qualities: vec![40; 50],
            },
            UniqueRecord {
                sequence: b"AAAAACCCCCGGGGGAAAAACACACGGGGGAAAAACCCCCGGGGGAAAAA".to_vec(),
                abundance: 200, // Very high abundance, should bud under loose omega_a
                average_qualities: vec![40.0; 50],
                consensus_qualities: vec![40; 50],
            },
            UniqueRecord {
                sequence: b"AAAAACCCCCGGGGGTATATCCCCCGGGGGAAAAACCCCCGGGGGAAAAA".to_vec(),
                abundance: 10,
                average_qualities: vec![40.0; 50],
                consensus_qualities: vec![40; 50],
            },
        ];

        let err_mat = generate_default_error_matrix(40);

        // Run with extremely stringent threshold (omega_a = 0.0) - should NOT bud, stay as 1 cluster
        let clustering_strict = Clustering::run_dada(
            &uniques,
            &err_mat,
            0.0,
            0.0,
            true,
            5,
            -4,
            -8,
            -2,
            16,
            false,
            false,
        );

        // Under extremely strict threshold, everything stays in 1 cluster
        assert_eq!(clustering_strict.clusters.len(), 1);

        // Run with loose threshold (omega_a = 1e-1) - should bud Rank 2 and Rank 3!
        println!("=== DEBUG: RUNNING LOOSE CLUSTERING ===");
        let clustering_loose = Clustering::run_dada(
            &uniques,
            &err_mat,
            1e-1,
            1e-1,
            true,
            5,
            -4,
            -8,
            -2,
            16,
            false,
            true, // Enable verbose
        );

        for (idx, record) in clustering_loose.records.iter().enumerate() {
            println!(
                "Record {}: abundance={}, p_value={:.2e}, lock={}, comp={:?}, e_minmax={:.4}",
                idx, record.abundance, record.p_value, record.lock, record.comp, record.e_minmax
            );
        }

        // Under loose threshold, it should successfully bud new clusters
        assert!(clustering_loose.clusters.len() > 1);
    }
}
