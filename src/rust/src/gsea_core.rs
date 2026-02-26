//! Core GSEA Enrichment Score calculation
//!
//! Implements the weighted Kolmogorov-Smirnov statistic for GSEA.
//!
//! Two implementations available:
//! - `calc_enrichment_score`: O(n) standard implementation
//! - `calc_enrichment_score_sparse`: O(k log k) sparse implementation (faster when k << n)

/// Result structure for a single pathway GSEA analysis
#[derive(Debug, Clone)]
pub struct GseaResult {
    pub pathway: String,
    pub es: f64,
    pub nes: f64,
    pub pval: f64,
    pub size: usize,
    pub leading_edge: Vec<usize>,
}

/// Calculate Enrichment Score using weighted Kolmogorov-Smirnov statistic
///
/// # Arguments
/// * `stats` - Gene statistics, should be sorted in decreasing order
/// * `gene_set_indices` - Indices of genes belonging to the gene set (0-based)
/// * `gsea_param` - Weighting exponent for gene statistics (typically 1.0)
///
/// # Returns
/// Tuple of (ES, leading_edge_indices)
pub fn calc_enrichment_score(
    stats: &[f64],
    gene_set_indices: &[usize],
    gsea_param: f64,
) -> (f64, Vec<usize>) {
    let n = stats.len();
    let k = gene_set_indices.len();

    if k == 0 || n == 0 {
        return (0.0, Vec::new());
    }

    // Create a set for O(1) lookup
    let gene_set: std::collections::HashSet<usize> =
        gene_set_indices.iter().cloned().collect();

    // Calculate the sum of weighted statistics for genes in the set
    let mut weight_sum = 0.0;
    for &idx in gene_set_indices {
        if idx < n {
            weight_sum += stats[idx].abs().powf(gsea_param);
        }
    }

    if weight_sum == 0.0 {
        return (0.0, Vec::new());
    }

    // Penalty for genes not in the set
    let miss_penalty = 1.0 / (n - k) as f64;

    // Running sum calculation
    let mut running_sum = 0.0;
    let mut max_es = 0.0;
    let mut min_es = 0.0;
    let mut max_pos = 0;
    let mut min_pos = 0;

    // Track leading edge
    let mut leading_edge_positive: Vec<usize> = Vec::new();
    let mut leading_edge_negative: Vec<usize> = Vec::new();

    for i in 0..n {
        if gene_set.contains(&i) {
            // Gene is in the set: add weighted contribution
            let weight = stats[i].abs().powf(gsea_param) / weight_sum;
            running_sum += weight;

            // Track genes contributing to positive enrichment
            if running_sum > 0.0 && running_sum > max_es {
                leading_edge_positive.push(i);
            }
            // Track genes contributing to negative enrichment
            if running_sum < 0.0 && running_sum < min_es {
                leading_edge_negative.push(i);
            }
        } else {
            // Gene is not in the set: subtract penalty
            running_sum -= miss_penalty;
        }

        // Track max and min
        if running_sum > max_es {
            max_es = running_sum;
            max_pos = i;
        }
        if running_sum < min_es {
            min_es = running_sum;
            min_pos = i;
        }
    }

    // ES is the maximum absolute deviation
    let (es, leading_edge) = if max_es.abs() >= min_es.abs() {
        // Positive enrichment: genes at top of list
        let le: Vec<usize> = gene_set_indices
            .iter()
            .filter(|&&idx| idx <= max_pos)
            .cloned()
            .collect();
        (max_es, le)
    } else {
        // Negative enrichment: genes at bottom of list
        let le: Vec<usize> = gene_set_indices
            .iter()
            .filter(|&&idx| idx >= min_pos)
            .cloned()
            .collect();
        (min_es, le)
    };

    (es, leading_edge)
}

/// Calculate ES for a random gene set of given size
/// Used for null distribution estimation
pub fn calc_random_es(stats: &[f64], set_size: usize, indices: &[usize], gsea_param: f64) -> f64 {
    let (es, _) = calc_enrichment_score_sparse(stats, &indices[..set_size], gsea_param);
    es
}

/// Batch calculate ES for multiple random gene sets
/// Optimized for Gamma fitting - calculates ES only (no leading edge)
///
/// # Arguments
/// * `stats` - Gene statistics (sorted)
/// * `set_size` - Size of each gene set
/// * `samples` - Vector of pre-sampled index arrays, each of length >= set_size
/// * `gsea_param` - GSEA weighting parameter
///
/// # Returns
/// Vector of ES values
#[inline]
pub fn calc_batch_es(
    stats: &[f64],
    set_size: usize,
    samples: &[Vec<usize>],
    gsea_param: f64,
) -> Vec<f64> {
    let n = stats.len();
    if set_size == 0 || n == 0 || samples.is_empty() {
        return vec![0.0; samples.len()];
    }

    let n_miss = n - set_size;
    let miss_penalty = if n_miss > 0 { 1.0 / n_miss as f64 } else { 0.0 };

    // Pre-compute absolute values raised to power (common for all samples)
    // This is a key optimization - avoid repeated powf calls
    let abs_stats_pow: Vec<f64> = if gsea_param == 1.0 {
        stats.iter().map(|x| x.abs()).collect()
    } else {
        stats.iter().map(|x| x.abs().powf(gsea_param)).collect()
    };

    samples
        .iter()
        .map(|sample| {
            calc_es_fast(&abs_stats_pow, set_size, sample, miss_penalty)
        })
        .collect()
}

/// Fast ES calculation using pre-computed absolute powers
/// Optimized inner loop for batch processing
#[inline(always)]
fn calc_es_fast(
    abs_stats_pow: &[f64],
    set_size: usize,
    indices: &[usize],
    miss_penalty: f64,
) -> f64 {
    let n = abs_stats_pow.len();

    // Sort indices for sparse traversal
    let mut sorted: Vec<usize> = indices[..set_size]
        .iter()
        .filter(|&&i| i < n)
        .cloned()
        .collect();
    sorted.sort_unstable();

    if sorted.is_empty() {
        return 0.0;
    }

    // Calculate weight sum
    let weight_sum: f64 = sorted.iter().map(|&i| abs_stats_pow[i]).sum();
    if weight_sum == 0.0 {
        return 0.0;
    }

    // Sparse ES calculation
    let mut running_sum = 0.0;
    let mut max_es = 0.0f64;
    let mut min_es = 0.0f64;
    let mut prev_pos = 0usize;

    for &idx in &sorted {
        // Accumulate misses
        let n_miss = idx.saturating_sub(prev_pos);
        if n_miss > 0 {
            running_sum -= miss_penalty * n_miss as f64;
            min_es = min_es.min(running_sum);
        }

        // Add hit
        running_sum += abs_stats_pow[idx] / weight_sum;
        max_es = max_es.max(running_sum);

        prev_pos = idx + 1;
    }

    // Trailing misses
    if prev_pos < n {
        running_sum -= miss_penalty * (n - prev_pos) as f64;
        min_es = min_es.min(running_sum);
    }

    // Return ES with larger absolute value
    if max_es.abs() >= min_es.abs() {
        max_es
    } else {
        min_es
    }
}

/// Sparse Enrichment Score calculation - O(k log k) instead of O(n)
///
/// This implementation is much faster when the gene set size k is much smaller
/// than the total number of genes n (which is the typical case).
///
/// # Algorithm
/// Instead of iterating through all n genes, we:
/// 1. Sort the gene set indices
/// 2. Only visit positions where genes are in the set
/// 3. Calculate cumulative miss penalty between hits using arithmetic
///
/// # Key Insight
/// Between two consecutive hits, running_sum only decreases (all misses).
/// Therefore:
/// - Maximum ES occurs immediately after a hit (just added weight)
/// - Minimum ES occurs immediately before a hit (after consecutive misses)
///
/// # Arguments
/// * `stats` - Gene statistics, should be sorted in decreasing order
/// * `gene_set_indices` - Indices of genes belonging to the gene set (0-based)
/// * `gsea_param` - Weighting exponent for gene statistics (typically 1.0)
///
/// # Returns
/// Tuple of (ES, leading_edge_indices)
pub fn calc_enrichment_score_sparse(
    stats: &[f64],
    gene_set_indices: &[usize],
    gsea_param: f64,
) -> (f64, Vec<usize>) {
    let n = stats.len();
    let k = gene_set_indices.len();

    if k == 0 || n == 0 {
        return (0.0, Vec::new());
    }

    // Sort gene set indices for sparse traversal
    let mut sorted_indices: Vec<usize> = gene_set_indices.to_vec();
    sorted_indices.sort_unstable();

    // Filter out invalid indices and calculate weight sum
    let mut weight_sum = 0.0;
    let mut valid_indices: Vec<usize> = Vec::with_capacity(k);

    for &idx in &sorted_indices {
        if idx < n {
            weight_sum += stats[idx].abs().powf(gsea_param);
            valid_indices.push(idx);
        }
    }

    let k_valid = valid_indices.len();
    if k_valid == 0 || weight_sum == 0.0 {
        return (0.0, Vec::new());
    }

    // Miss penalty per non-hit gene
    let miss_penalty = 1.0 / (n - k_valid) as f64;

    // Sparse running sum calculation
    let mut running_sum = 0.0;
    let mut max_es = 0.0;
    let mut min_es = 0.0;
    let mut max_idx = 0usize;  // Index where max ES occurred
    let mut min_idx = 0usize;  // Index where min ES occurred
    let mut prev_pos = 0usize; // Position after previous hit (or 0 at start)

    for &idx in &valid_indices {
        // Calculate misses between previous position and current hit
        let n_miss = idx.saturating_sub(prev_pos);

        if n_miss > 0 {
            // Accumulate miss penalty
            running_sum -= miss_penalty * n_miss as f64;

            // Check for new minimum (occurs just before this hit)
            if running_sum < min_es {
                min_es = running_sum;
                min_idx = idx.saturating_sub(1);
            }
        }

        // Add hit contribution
        let weight = stats[idx].abs().powf(gsea_param) / weight_sum;
        running_sum += weight;

        // Check for new maximum (occurs just after this hit)
        if running_sum > max_es {
            max_es = running_sum;
            max_idx = idx;
        }

        prev_pos = idx + 1;
    }

    // Handle trailing misses after last hit
    if prev_pos < n {
        let n_miss = n - prev_pos;
        running_sum -= miss_penalty * n_miss as f64;

        if running_sum < min_es {
            min_es = running_sum;
            min_idx = n - 1;
        }
    }

    // Determine ES and leading edge
    let (es, leading_edge) = if max_es.abs() >= min_es.abs() {
        // Positive enrichment: genes at top of list (up to max position)
        let le: Vec<usize> = valid_indices
            .iter()
            .filter(|&&idx| idx <= max_idx)
            .cloned()
            .collect();
        (max_es, le)
    } else {
        // Negative enrichment: genes at bottom of list (from min position)
        let le: Vec<usize> = valid_indices
            .iter()
            .filter(|&&idx| idx >= min_idx)
            .cloned()
            .collect();
        (min_es, le)
    };

    (es, leading_edge)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_gene_set() {
        let stats = vec![1.0, 0.5, 0.0, -0.5, -1.0];
        let indices: Vec<usize> = vec![];
        let (es, le) = calc_enrichment_score(&stats, &indices, 1.0);
        assert_eq!(es, 0.0);
        assert!(le.is_empty());

        // Sparse version should match
        let (es_sparse, le_sparse) = calc_enrichment_score_sparse(&stats, &indices, 1.0);
        assert_eq!(es_sparse, 0.0);
        assert!(le_sparse.is_empty());
    }

    #[test]
    fn test_top_enriched() {
        // Gene set at the top of the ranked list should give positive ES
        let stats = vec![2.0, 1.5, 1.0, 0.5, 0.0, -0.5, -1.0, -1.5, -2.0];
        let indices = vec![0, 1, 2]; // Top 3 genes
        let (es, _) = calc_enrichment_score(&stats, &indices, 1.0);
        assert!(es > 0.0, "Expected positive ES, got {}", es);

        // Sparse version should match
        let (es_sparse, _) = calc_enrichment_score_sparse(&stats, &indices, 1.0);
        assert!((es - es_sparse).abs() < 1e-10, "Sparse ES {} != Original ES {}", es_sparse, es);
    }

    #[test]
    fn test_bottom_enriched() {
        // Gene set at the bottom of the ranked list should give negative ES
        let stats = vec![2.0, 1.5, 1.0, 0.5, 0.0, -0.5, -1.0, -1.5, -2.0];
        let indices = vec![6, 7, 8]; // Bottom 3 genes
        let (es, _) = calc_enrichment_score(&stats, &indices, 1.0);
        assert!(es < 0.0, "Expected negative ES, got {}", es);

        // Sparse version should match
        let (es_sparse, _) = calc_enrichment_score_sparse(&stats, &indices, 1.0);
        assert!((es - es_sparse).abs() < 1e-10, "Sparse ES {} != Original ES {}", es_sparse, es);
    }

    #[test]
    fn test_sparse_vs_original_consistency() {
        // Test with various gene set configurations
        let n = 1000;
        let stats: Vec<f64> = (0..n).map(|i| 3.0 - (i as f64 / 166.0)).collect();

        // Test cases: different sizes and positions
        let test_cases = vec![
            (0..50).collect::<Vec<_>>(),           // Top 50
            (950..1000).collect::<Vec<_>>(),       // Bottom 50
            (475..525).collect::<Vec<_>>(),        // Middle 50
            vec![0, 100, 200, 300, 400, 500, 600, 700, 800, 900], // Spread out
            vec![0, 1, 998, 999],                  // Extremes
            (0..200).collect::<Vec<_>>(),          // Larger set
        ];

        for indices in test_cases {
            let (es_orig, le_orig) = calc_enrichment_score(&stats, &indices, 1.0);
            let (es_sparse, le_sparse) = calc_enrichment_score_sparse(&stats, &indices, 1.0);

            assert!(
                (es_orig - es_sparse).abs() < 1e-10,
                "ES mismatch for {:?}: original={}, sparse={}",
                &indices[..indices.len().min(5)],
                es_orig,
                es_sparse
            );

            // Leading edge should have same genes (order may differ)
            let mut le_orig_sorted = le_orig.clone();
            let mut le_sparse_sorted = le_sparse.clone();
            le_orig_sorted.sort();
            le_sparse_sorted.sort();
            assert_eq!(
                le_orig_sorted, le_sparse_sorted,
                "Leading edge mismatch for indices starting with {:?}",
                &indices[..indices.len().min(5)]
            );
        }
    }

    #[test]
    fn test_sparse_with_unsorted_indices() {
        // Sparse version should handle unsorted input
        let stats = vec![2.0, 1.5, 1.0, 0.5, 0.0, -0.5, -1.0, -1.5, -2.0];
        let indices_unsorted = vec![2, 0, 1]; // Unsorted
        let indices_sorted = vec![0, 1, 2];   // Sorted

        let (es_unsorted, _) = calc_enrichment_score_sparse(&stats, &indices_unsorted, 1.0);
        let (es_sorted, _) = calc_enrichment_score_sparse(&stats, &indices_sorted, 1.0);

        assert!((es_unsorted - es_sorted).abs() < 1e-10);
    }

    #[test]
    fn test_sparse_performance_scenario() {
        // Simulate typical GSEA scenario: n=20000 genes, k=50 gene set
        let n = 20000;
        let k = 50;
        let stats: Vec<f64> = (0..n).map(|i| 5.0 - (i as f64 / 2000.0)).collect();
        let indices: Vec<usize> = (0..k).map(|i| i * (n / k)).collect();

        let (es_orig, _) = calc_enrichment_score(&stats, &indices, 1.0);
        let (es_sparse, _) = calc_enrichment_score_sparse(&stats, &indices, 1.0);

        assert!(
            (es_orig - es_sparse).abs() < 1e-10,
            "Large scale test failed: original={}, sparse={}",
            es_orig,
            es_sparse
        );
    }
}
