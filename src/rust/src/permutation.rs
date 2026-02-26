//! Fallback permutation-based p-value calculation
//!
//! Used when Gamma distribution fitting fails (e.g., very small gene sets).

use rand::seq::SliceRandom;
use rand::thread_rng;

use crate::gsea_core::calc_enrichment_score_sparse;

/// Calculate p-value using permutation testing
///
/// This is the fallback method when Gamma fitting fails.
/// Generates random gene sets and calculates empirical p-value.
///
/// # Arguments
/// * `stats` - Gene statistics (sorted)
/// * `gene_set_indices` - Indices of genes in the target set
/// * `gsea_param` - GSEA weighting parameter
/// * `n_perm` - Number of permutations
///
/// # Returns
/// Tuple of (p-value, mean_es_positive, mean_es_negative) for NES calculation
pub fn permutation_pvalue(
    stats: &[f64],
    gene_set_indices: &[usize],
    gsea_param: f64,
    n_perm: usize,
) -> (f64, f64, f64) {
    let n = stats.len();
    let k = gene_set_indices.len();

    if k == 0 || n == 0 {
        return (1.0, 0.0, 0.0);
    }

    // Calculate observed ES
    let (observed_es, _) = calc_enrichment_score_sparse(stats, gene_set_indices, gsea_param);

    if observed_es == 0.0 {
        return (1.0, 0.0, 0.0);
    }

    let mut rng = thread_rng();
    let all_indices: Vec<usize> = (0..n).collect();

    let mut more_extreme = 0;
    let mut positive_es_sum = 0.0;
    let mut negative_es_sum = 0.0;
    let mut positive_count = 0;
    let mut negative_count = 0;

    for _ in 0..n_perm {
        // Random sample of gene indices
        let sample: Vec<usize> = all_indices
            .choose_multiple(&mut rng, k)
            .cloned()
            .collect();

        let (random_es, _) = calc_enrichment_score_sparse(stats, &sample, gsea_param);

        // Track for NES calculation
        if random_es > 0.0 {
            positive_es_sum += random_es;
            positive_count += 1;
        } else if random_es < 0.0 {
            negative_es_sum += random_es.abs();
            negative_count += 1;
        }

        // Count more extreme values
        if observed_es > 0.0 {
            if random_es >= observed_es {
                more_extreme += 1;
            }
        } else {
            if random_es <= observed_es {
                more_extreme += 1;
            }
        }
    }

    // Calculate p-value with pseudocount to avoid exact 0
    let pval = (more_extreme as f64 + 1.0) / (n_perm as f64 + 1.0);

    // Calculate mean absolute ES for normalization
    let mean_pos = if positive_count > 0 {
        positive_es_sum / positive_count as f64
    } else {
        0.0
    };

    let mean_neg = if negative_count > 0 {
        negative_es_sum / negative_count as f64
    } else {
        0.0
    };

    (pval, mean_pos, mean_neg)
}

/// Calculate NES from permutation results
pub fn calculate_nes_permutation(es: f64, mean_pos: f64, mean_neg: f64) -> f64 {
    if es == 0.0 {
        return 0.0;
    }

    if es > 0.0 {
        if mean_pos > 0.0 {
            es / mean_pos
        } else {
            es
        }
    } else {
        if mean_neg > 0.0 {
            es / mean_neg // Note: es is negative, mean_neg is positive
        } else {
            es
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_permutation_pvalue() {
        // Create sorted stats
        let stats: Vec<f64> = (0..100)
            .map(|i| 1.0 - (i as f64 / 50.0))
            .collect();

        // Gene set at top should have low p-value
        let top_genes: Vec<usize> = (0..10).collect();
        let (pval_top, _, _) = permutation_pvalue(&stats, &top_genes, 1.0, 1000);

        // Gene set in middle should have high p-value
        let mid_genes: Vec<usize> = (45..55).collect();
        let (pval_mid, _, _) = permutation_pvalue(&stats, &mid_genes, 1.0, 1000);

        assert!(pval_top < pval_mid, "Top genes p-value {} should be < middle genes p-value {}", pval_top, pval_mid);
    }

    #[test]
    fn test_nes_calculation() {
        let nes = calculate_nes_permutation(0.5, 0.2, 0.2);
        assert!((nes - 2.5).abs() < 0.01);

        let nes_neg = calculate_nes_permutation(-0.5, 0.2, 0.2);
        assert!((nes_neg - (-2.5)).abs() < 0.01);
    }
}
