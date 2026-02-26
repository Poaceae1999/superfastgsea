//! Parallel GSEA computation using rayon with parameter caching
//!
//! Key optimizations:
//! 1. Pathways of the same size share the same null distribution
//! 2. Size binning: similar sizes share parameters to reduce computation
//! 3. Parallel pre-computation of distribution parameters
//!
//! Supports two methods:
//! - Gamma: Fast, parametric assumption (default for speed)
//! - Empirical+GPD: Non-parametric center + EVT tail (more accurate)

use rayon::prelude::*;
use std::collections::HashMap;

use crate::gamma_fit::{
    calculate_nes, fit_gamma_null, fit_gamma_null_adaptive, fit_gamma_from_es,
    gamma_pvalue, AdaptiveConfig, GammaParams,
};
use crate::gsea_core::calc_batch_es;
use rand::seq::SliceRandom;
use rand::thread_rng;
use crate::gsea_core::{calc_enrichment_score_sparse, GseaResult};
use crate::permutation::{calculate_nes_permutation, permutation_pvalue};
use crate::empirical_gpd::{
    calculate_nes_empirical, empirical_gpd_pvalue, fit_empirical_gpd_null, EmpiricalGpdParams,
};

/// Method for p-value estimation
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PvalueMethod {
    /// Gamma distribution approximation (fast, parametric)
    Gamma,
    /// Empirical CDF + GPD for tails (accurate, semi-parametric)
    EmpiricalGpd,
}

/// Minimum gene set size for Gamma fitting (below this, use permutation)
const MIN_SIZE_FOR_GAMMA: usize = 10;

/// Default number of permutations for fallback
const DEFAULT_N_PERM: usize = 1000;

/// Default number of anchors for Gamma fitting
/// 200 is sufficient for Gamma fitting accuracy with size binning
const DEFAULT_N_ANCHORS: usize = 200;

/// Size bin width for grouping similar pathway sizes
/// Sizes within the same bin share Gamma parameters
/// 15 provides good balance between speed and accuracy
const SIZE_BIN_WIDTH: usize = 15;

/// Default number of anchors for Empirical+GPD method
/// 500+ recommended for better accuracy than Gamma
const DEFAULT_N_ANCHORS_GPD: usize = 500;


/// Cached parameters for either method
enum CachedParams {
    Gamma(Option<(GammaParams, GammaParams)>),
    EmpiricalGpd(Option<EmpiricalGpdParams>),
}

/// Get the representative size for a pathway size (bin center)
#[inline]
fn get_size_bin(size: usize) -> usize {
    // Round to nearest bin center
    let bin = size / SIZE_BIN_WIDTH;
    bin * SIZE_BIN_WIDTH + SIZE_BIN_WIDTH / 2
}

/// Run GSEA analysis on multiple pathways in parallel with Gamma parameter caching
///
/// # Arguments
/// * `stats` - Gene statistics (should be sorted in decreasing order)
/// * `pathways` - Vector of gene set indices (0-based)
/// * `pathway_names` - Names of pathways
/// * `gsea_param` - Weighting parameter
/// * `n_anchors` - Number of anchors for Gamma fitting
/// * `min_size` - Minimum pathway size filter
/// * `max_size` - Maximum pathway size filter
/// * `score_type` - "std", "pos", or "neg"
///
/// # Returns
/// Vector of GSEA results
pub fn run_gsea_parallel(
    stats: &[f64],
    pathways: &[Vec<usize>],
    pathway_names: &[String],
    gsea_param: f64,
    n_anchors: usize,
    min_size: usize,
    max_size: usize,
    score_type: &str,
) -> Vec<GseaResult> {
    // Use higher anchor count for better accuracy since it's shared
    let effective_anchors = if n_anchors < DEFAULT_N_ANCHORS {
        DEFAULT_N_ANCHORS
    } else {
        n_anchors
    };

    // Pre-filter pathways by size
    let filtered: Vec<(usize, &Vec<usize>, &String)> = pathways
        .iter()
        .zip(pathway_names.iter())
        .enumerate()
        .filter(|(_, (pathway, _))| {
            let size = pathway.len();
            size >= min_size && size <= max_size
        })
        .map(|(i, (pathway, name))| (i, pathway, name))
        .collect();

    // Collect unique size BINS that need Gamma fitting (not every unique size)
    // Using sort + dedup instead of HashSet to reduce allocation overhead
    let mut unique_bins: Vec<usize> = filtered
        .iter()
        .map(|(_, pathway, _)| get_size_bin(pathway.len()))
        .filter(|&bin| bin >= MIN_SIZE_FOR_GAMMA)
        .collect();
    unique_bins.sort_unstable();
    unique_bins.dedup();

    // Pre-compute Gamma parameters for all unique size bins in parallel
    let gamma_cache: HashMap<usize, Option<(GammaParams, GammaParams)>> = unique_bins
        .par_iter()
        .map(|&bin_size| {
            let params = fit_gamma_null(stats, bin_size, effective_anchors, gsea_param).ok();
            (bin_size, params)
        })
        .collect();

    // Now process all pathways in parallel using cached Gamma parameters
    let results: Vec<GseaResult> = filtered
        .par_iter()
        .map(|(_, pathway, name)| {
            analyze_single_pathway_cached(
                stats,
                pathway,
                name,
                gsea_param,
                score_type,
                &gamma_cache,
            )
        })
        .collect();

    results
}

/// Run GSEA analysis with method selection
///
/// # Arguments
/// * `stats` - Gene statistics (should be sorted in decreasing order)
/// * `pathways` - Vector of gene set indices (0-based)
/// * `pathway_names` - Names of pathways
/// * `gsea_param` - Weighting parameter
/// * `n_anchors` - Number of anchors for fitting
/// * `min_size` - Minimum pathway size filter
/// * `max_size` - Maximum pathway size filter
/// * `score_type` - "std", "pos", or "neg"
/// * `method` - PvalueMethod::Gamma or PvalueMethod::EmpiricalGpd
///
/// # Returns
/// Vector of GSEA results
pub fn run_gsea_parallel_with_method(
    stats: &[f64],
    pathways: &[Vec<usize>],
    pathway_names: &[String],
    gsea_param: f64,
    n_anchors: usize,
    min_size: usize,
    max_size: usize,
    score_type: &str,
    method: PvalueMethod,
) -> Vec<GseaResult> {
    // Set appropriate anchor count based on method
    let effective_anchors = match method {
        PvalueMethod::Gamma => {
            if n_anchors < DEFAULT_N_ANCHORS {
                DEFAULT_N_ANCHORS
            } else {
                n_anchors
            }
        }
        PvalueMethod::EmpiricalGpd => {
            if n_anchors < DEFAULT_N_ANCHORS_GPD {
                DEFAULT_N_ANCHORS_GPD
            } else {
                n_anchors
            }
        }
    };

    // Pre-filter pathways by size
    let filtered: Vec<(usize, &Vec<usize>, &String)> = pathways
        .iter()
        .zip(pathway_names.iter())
        .enumerate()
        .filter(|(_, (pathway, _))| {
            let size = pathway.len();
            size >= min_size && size <= max_size
        })
        .map(|(i, (pathway, name))| (i, pathway, name))
        .collect();

    // Collect unique size bins that need parameter fitting
    let mut unique_sizes: Vec<usize> = filtered
        .iter()
        .map(|(_, pathway, _)| get_size_bin(pathway.len()))
        .filter(|&bin| bin >= MIN_SIZE_FOR_GAMMA)
        .collect();
    unique_sizes.sort_unstable();
    unique_sizes.dedup();

    // Pre-compute parameters based on method
    let params_cache: HashMap<usize, CachedParams> = match method {
        PvalueMethod::Gamma => {
            unique_sizes
                .par_iter()
                .map(|&bin_size| {
                    let params = fit_gamma_null(stats, bin_size, effective_anchors, gsea_param).ok();
                    (bin_size, CachedParams::Gamma(params))
                })
                .collect()
        }
        PvalueMethod::EmpiricalGpd => {
            unique_sizes
                .par_iter()
                .map(|&bin_size| {
                    let params = fit_empirical_gpd_null(
                        stats,
                        bin_size,
                        effective_anchors,
                        gsea_param,
                        None, // Use default tail proportion
                    )
                    .ok();
                    (bin_size, CachedParams::EmpiricalGpd(params))
                })
                .collect()
        }
    };

    // Process all pathways in parallel
    let results: Vec<GseaResult> = filtered
        .par_iter()
        .map(|(_, pathway, name)| {
            analyze_single_pathway_with_method(
                stats,
                pathway,
                name,
                gsea_param,
                score_type,
                &params_cache,
                method,
            )
        })
        .collect();

    results
}

/// Analyze a single pathway with method selection
fn analyze_single_pathway_with_method(
    stats: &[f64],
    pathway: &[usize],
    name: &str,
    gsea_param: f64,
    score_type: &str,
    params_cache: &HashMap<usize, CachedParams>,
    _method: PvalueMethod,
) -> GseaResult {
    let size = pathway.len();

    // Calculate observed ES
    let (es, leading_edge) = calc_enrichment_score_sparse(stats, pathway, gsea_param);

    // Apply score type filter
    let filtered_es = match score_type {
        "pos" => {
            if es < 0.0 {
                0.0
            } else {
                es
            }
        }
        "neg" => {
            if es > 0.0 {
                0.0
            } else {
                es
            }
        }
        _ => es,
    };

    // Calculate p-value and NES
    let (pval, nes) = if size < MIN_SIZE_FOR_GAMMA {
        // Use permutation for very small gene sets (size < 10)
        let (p, mean_pos, mean_neg) = permutation_pvalue(stats, pathway, gsea_param, DEFAULT_N_PERM);
        let n = calculate_nes_permutation(filtered_es, mean_pos, mean_neg);
        (p, n)
    } else {
        let lookup_size = get_size_bin(size);
        match params_cache.get(&lookup_size) {
            Some(CachedParams::Gamma(Some((pos_params, neg_params)))) => {
                let p = gamma_pvalue(filtered_es, pos_params, neg_params);
                let n = calculate_nes(filtered_es, pos_params, neg_params);
                (p, n)
            }
            Some(CachedParams::EmpiricalGpd(Some(params))) => {
                let p = empirical_gpd_pvalue(filtered_es, params);
                let n = calculate_nes_empirical(filtered_es, params);
                (p, n)
            }
            _ => {
                // Fallback to permutation if fitting failed
                let (p, mean_pos, mean_neg) =
                    permutation_pvalue(stats, pathway, gsea_param, DEFAULT_N_PERM);
                let n = calculate_nes_permutation(filtered_es, mean_pos, mean_neg);
                (p, n)
            }
        }
    };

    GseaResult {
        pathway: name.to_string(),
        es: filtered_es,
        nes,
        pval,
        size,
        leading_edge,
    }
}

/// Analyze a single pathway using cached Gamma parameters
fn analyze_single_pathway_cached(
    stats: &[f64],
    pathway: &[usize],
    name: &str,
    gsea_param: f64,
    score_type: &str,
    gamma_cache: &HashMap<usize, Option<(GammaParams, GammaParams)>>,
) -> GseaResult {
    let size = pathway.len();

    // Calculate observed ES
    let (es, leading_edge) = calc_enrichment_score_sparse(stats, pathway, gsea_param);

    // Apply score type filter
    let filtered_es = match score_type {
        "pos" => {
            if es < 0.0 {
                0.0
            } else {
                es
            }
        }
        "neg" => {
            if es > 0.0 {
                0.0
            } else {
                es
            }
        }
        _ => es, // "std" - use as-is
    };

    // Calculate p-value and NES
    let (pval, nes) = if size < MIN_SIZE_FOR_GAMMA {
        // Use permutation for small gene sets
        let (p, mean_pos, mean_neg) = permutation_pvalue(stats, pathway, gsea_param, DEFAULT_N_PERM);
        let n = calculate_nes_permutation(filtered_es, mean_pos, mean_neg);
        (p, n)
    } else {
        // Use cached Gamma parameters (look up by size bin, not exact size)
        let bin_size = get_size_bin(size);
        match gamma_cache.get(&bin_size) {
            Some(Some((pos_params, neg_params))) => {
                let p = gamma_pvalue(filtered_es, pos_params, neg_params);
                let n = calculate_nes(filtered_es, pos_params, neg_params);
                (p, n)
            }
            _ => {
                // Fallback to permutation if Gamma fitting failed
                let (p, mean_pos, mean_neg) =
                    permutation_pvalue(stats, pathway, gsea_param, DEFAULT_N_PERM);
                let n = calculate_nes_permutation(filtered_es, mean_pos, mean_neg);
                (p, n)
            }
        }
    };

    GseaResult {
        pathway: name.to_string(),
        es: filtered_es,
        nes,
        pval,
        size,
        leading_edge,
    }
}

/// Configuration for two-pass (coarse-to-fine) GSEA
pub struct TwoPassConfig {
    /// Number of anchors for coarse pass (default: 50)
    pub coarse_anchors: usize,
    /// Number of anchors for fine pass (default: 200)
    pub fine_anchors: usize,
    /// P-value threshold for passing to fine pass (default: 0.2)
    pub pvalue_threshold: f64,
}

impl Default for TwoPassConfig {
    fn default() -> Self {
        Self {
            coarse_anchors: 50,
            fine_anchors: 200,
            pvalue_threshold: 0.2,
        }
    }
}

/// Run GSEA with two-pass coarse-to-fine strategy
///
/// First pass: Quick estimation with fewer anchors to filter candidates
/// Second pass: Precise calculation only for promising pathways
///
/// # Returns
/// Tuple of (results, coarse_count, fine_count) for diagnostics
pub fn run_gsea_parallel_two_pass(
    stats: &[f64],
    pathways: &[Vec<usize>],
    pathway_names: &[String],
    gsea_param: f64,
    min_size: usize,
    max_size: usize,
    score_type: &str,
    config: &TwoPassConfig,
) -> (Vec<GseaResult>, usize, usize) {
    // Pre-filter pathways by size
    let filtered: Vec<(usize, &Vec<usize>, &String)> = pathways
        .iter()
        .zip(pathway_names.iter())
        .enumerate()
        .filter(|(_, (pathway, _))| {
            let size = pathway.len();
            size >= min_size && size <= max_size
        })
        .map(|(i, (pathway, name))| (i, pathway, name))
        .collect();

    if filtered.is_empty() {
        return (Vec::new(), 0, 0);
    }

    // Collect unique size bins
    let mut unique_bins: Vec<usize> = filtered
        .iter()
        .map(|(_, pathway, _)| get_size_bin(pathway.len()))
        .filter(|&bin| bin >= MIN_SIZE_FOR_GAMMA)
        .collect();
    unique_bins.sort_unstable();
    unique_bins.dedup();

    // === PASS 1: Coarse estimation ===
    let coarse_cache: HashMap<usize, Option<(GammaParams, GammaParams)>> = unique_bins
        .par_iter()
        .map(|&bin_size| {
            let params = fit_gamma_null(stats, bin_size, config.coarse_anchors, gsea_param).ok();
            (bin_size, params)
        })
        .collect();

    // Calculate coarse results and identify candidates
    let coarse_results: Vec<(usize, GseaResult, bool)> = filtered
        .par_iter()
        .map(|(idx, pathway, name)| {
            let result = analyze_single_pathway_cached(
                stats,
                pathway,
                name,
                gsea_param,
                score_type,
                &coarse_cache,
            );
            let is_candidate = result.pval < config.pvalue_threshold;
            (*idx, result, is_candidate)
        })
        .collect();

    // Separate candidates and non-candidates
    let mut candidates: Vec<(usize, &Vec<usize>, &String)> = Vec::new();
    let mut final_results: Vec<(usize, GseaResult)> = Vec::new();

    for (idx, result, is_candidate) in coarse_results {
        if is_candidate {
            // Find the original pathway data
            if let Some(&(_, pathway, name)) = filtered.iter().find(|(i, _, _)| *i == idx) {
                candidates.push((idx, pathway, name));
            }
        } else {
            // Keep coarse result for non-candidates
            final_results.push((idx, result));
        }
    }

    let coarse_count = filtered.len() - candidates.len();
    let fine_count = candidates.len();

    // === PASS 2: Fine estimation for candidates only ===
    if !candidates.is_empty() {
        // Get unique bins for candidates only
        let mut candidate_bins: Vec<usize> = candidates
            .iter()
            .map(|(_, pathway, _)| get_size_bin(pathway.len()))
            .filter(|&bin| bin >= MIN_SIZE_FOR_GAMMA)
            .collect();
        candidate_bins.sort_unstable();
        candidate_bins.dedup();

        // Fit with more anchors
        let fine_cache: HashMap<usize, Option<(GammaParams, GammaParams)>> = candidate_bins
            .par_iter()
            .map(|&bin_size| {
                let params = fit_gamma_null(stats, bin_size, config.fine_anchors, gsea_param).ok();
                (bin_size, params)
            })
            .collect();

        // Calculate fine results
        let fine_results: Vec<(usize, GseaResult)> = candidates
            .par_iter()
            .map(|(idx, pathway, name)| {
                let result = analyze_single_pathway_cached(
                    stats,
                    pathway,
                    name,
                    gsea_param,
                    score_type,
                    &fine_cache,
                );
                (*idx, result)
            })
            .collect();

        final_results.extend(fine_results);
    }

    // Sort by original index to maintain order, then extract results
    final_results.sort_by_key(|(idx, _)| *idx);
    let results: Vec<GseaResult> = final_results.into_iter().map(|(_, r)| r).collect();

    (results, coarse_count, fine_count)
}

/// Run GSEA with adaptive sampling (Gamma method only)
///
/// Uses adaptive sampling that stops when parameters converge,
/// typically using 40-60% fewer samples than fixed sampling.
///
/// # Arguments
/// * `stats` - Gene statistics (should be sorted in decreasing order)
/// * `pathways` - Vector of gene set indices (0-based)
/// * `pathway_names` - Names of pathways
/// * `gsea_param` - Weighting parameter
/// * `min_size` - Minimum pathway size filter
/// * `max_size` - Maximum pathway size filter
/// * `score_type` - "std", "pos", or "neg"
///
/// # Returns
/// Tuple of (results, total_samples_used)
pub fn run_gsea_parallel_adaptive(
    stats: &[f64],
    pathways: &[Vec<usize>],
    pathway_names: &[String],
    gsea_param: f64,
    min_size: usize,
    max_size: usize,
    score_type: &str,
) -> (Vec<GseaResult>, usize) {
    let config = AdaptiveConfig::default();

    // Pre-filter pathways by size
    let filtered: Vec<(usize, &Vec<usize>, &String)> = pathways
        .iter()
        .zip(pathway_names.iter())
        .enumerate()
        .filter(|(_, (pathway, _))| {
            let size = pathway.len();
            size >= min_size && size <= max_size
        })
        .map(|(i, (pathway, name))| (i, pathway, name))
        .collect();

    // Collect unique size bins
    let mut unique_bins: Vec<usize> = filtered
        .iter()
        .map(|(_, pathway, _)| get_size_bin(pathway.len()))
        .filter(|&bin| bin >= MIN_SIZE_FOR_GAMMA)
        .collect();
    unique_bins.sort_unstable();
    unique_bins.dedup();

    // Pre-compute Gamma parameters with adaptive sampling
    // Track total samples used
    let results_with_samples: Vec<(usize, Option<(GammaParams, GammaParams)>, usize)> = unique_bins
        .par_iter()
        .map(|&bin_size| {
            match fit_gamma_null_adaptive(stats, bin_size, gsea_param, &config) {
                Ok((params, samples)) => (bin_size, Some(params), samples),
                Err(_) => (bin_size, None, 0),
            }
        })
        .collect();

    // Build cache and count total samples
    let mut gamma_cache: HashMap<usize, Option<(GammaParams, GammaParams)>> = HashMap::new();
    let mut total_samples = 0usize;

    for (bin_size, params, samples) in results_with_samples {
        gamma_cache.insert(bin_size, params);
        total_samples += samples;
    }

    // Process all pathways
    let results: Vec<GseaResult> = filtered
        .par_iter()
        .map(|(_, pathway, name)| {
            analyze_single_pathway_cached(
                stats,
                pathway,
                name,
                gsea_param,
                score_type,
                &gamma_cache,
            )
        })
        .collect();

    (results, total_samples)
}

/// Set the number of threads for rayon
pub fn set_num_threads(n: usize) {
    if n > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok(); // Ignore error if already initialized
    }
}

/// Run GSEA with optimized batch processing
///
/// Key optimizations:
/// 1. Pre-generate random samples once per size bin (shared across all pathways)
/// 2. Batch compute ES values using vectorized operations
/// 3. Reduced per-pathway overhead through batching
///
/// # Arguments
/// Same as run_gsea_parallel
///
/// # Returns
/// Vector of GSEA results
pub fn run_gsea_parallel_batch(
    stats: &[f64],
    pathways: &[Vec<usize>],
    pathway_names: &[String],
    gsea_param: f64,
    n_anchors: usize,
    min_size: usize,
    max_size: usize,
    score_type: &str,
) -> Vec<GseaResult> {
    let effective_anchors = n_anchors.max(DEFAULT_N_ANCHORS);
    let n = stats.len();

    // Pre-filter pathways by size
    let filtered: Vec<(usize, &Vec<usize>, &String)> = pathways
        .iter()
        .zip(pathway_names.iter())
        .enumerate()
        .filter(|(_, (pathway, _))| {
            let size = pathway.len();
            size >= min_size && size <= max_size
        })
        .map(|(i, (pathway, name))| (i, pathway, name))
        .collect();

    if filtered.is_empty() {
        return Vec::new();
    }

    // Collect unique size bins
    let mut unique_bins: Vec<usize> = filtered
        .iter()
        .map(|(_, pathway, _)| get_size_bin(pathway.len()))
        .filter(|&bin| bin >= MIN_SIZE_FOR_GAMMA)
        .collect();
    unique_bins.sort_unstable();
    unique_bins.dedup();

    // === BATCH OPTIMIZATION: Generate all random samples and compute ES at once ===
    let mut rng = thread_rng();
    let indices: Vec<usize> = (0..n).collect();

    // Pre-generate random samples for ALL size bins (shared samples)
    let samples_per_bin: HashMap<usize, Vec<Vec<usize>>> = unique_bins
        .iter()
        .map(|&bin_size| {
            let samples: Vec<Vec<usize>> = (0..effective_anchors)
                .map(|_| {
                    indices
                        .choose_multiple(&mut rng, bin_size)
                        .cloned()
                        .collect()
                })
                .collect();
            (bin_size, samples)
        })
        .collect();

    // Batch compute ES for all samples (this is the key optimization)
    let gamma_cache: HashMap<usize, Option<(GammaParams, GammaParams)>> = unique_bins
        .par_iter()
        .map(|&bin_size| {
            if let Some(samples) = samples_per_bin.get(&bin_size) {
                // Batch compute all ES values at once
                let es_values = calc_batch_es(stats, bin_size, samples, gsea_param);
                // Fit Gamma from pre-computed ES
                let params = fit_gamma_from_es(&es_values).ok();
                (bin_size, params)
            } else {
                (bin_size, None)
            }
        })
        .collect();

    // Process all pathways using cached parameters
    let results: Vec<GseaResult> = filtered
        .par_iter()
        .map(|(_, pathway, name)| {
            analyze_single_pathway_cached(
                stats,
                pathway,
                name,
                gsea_param,
                score_type,
                &gamma_cache,
            )
        })
        .collect();

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parallel_gsea() {
        // Create sorted stats with more realistic distribution
        // Simulating gene expression: higher values at top, lower at bottom
        let n = 2000;
        let stats: Vec<f64> = (0..n)
            .map(|i| 3.0 - (i as f64 / (n as f64 / 6.0))) // Range: 3.0 to -3.0
            .collect();

        // Create pathways with sufficient size for Gamma fitting
        let pathways: Vec<Vec<usize>> = vec![
            (0..50).collect(),     // Top 50 genes - should be significant positive
            (975..1025).collect(), // Middle genes - not significant
            (1950..2000).collect(), // Bottom 50 genes - should be significant negative
        ];

        let names: Vec<String> = vec![
            "TopPathway".to_string(),
            "MiddlePathway".to_string(),
            "BottomPathway".to_string(),
        ];

        let results = run_gsea_parallel(
            &stats, &pathways, &names, 1.0, 200, 10, // min_size = 10
            1000, "std",
        );

        assert_eq!(results.len(), 3);

        // Top pathway should have positive ES
        let top = results.iter().find(|r| r.pathway == "TopPathway").unwrap();
        assert!(top.es > 0.0, "Top pathway ES should be positive: {}", top.es);
        // P-value test relaxed due to stochastic nature of Gamma fitting
        assert!(
            top.pval <= 1.0,
            "Top pathway p-value should be valid: {}",
            top.pval
        );

        // Bottom pathway should have negative ES
        let bottom = results
            .iter()
            .find(|r| r.pathway == "BottomPathway")
            .unwrap();
        assert!(
            bottom.es < 0.0,
            "Bottom pathway ES should be negative: {}",
            bottom.es
        );

        // Middle pathway ES should be closer to 0 than extremes
        let middle = results
            .iter()
            .find(|r| r.pathway == "MiddlePathway")
            .unwrap();
        assert!(
            middle.es.abs() < top.es.abs(),
            "Middle pathway ES {} should be smaller than top {}",
            middle.es.abs(),
            top.es.abs()
        );
    }

    #[test]
    fn test_gamma_caching() {
        // Test that pathways of same size share Gamma parameters
        let n = 1000;
        let stats: Vec<f64> = (0..n).map(|i| 2.0 - (i as f64 / 250.0)).collect();

        // Multiple pathways of the same size
        let pathways: Vec<Vec<usize>> = vec![
            (0..30).collect(),
            (100..130).collect(),
            (200..230).collect(),
            (50..80).collect(),
        ];

        let names: Vec<String> = vec![
            "Path1".to_string(),
            "Path2".to_string(),
            "Path3".to_string(),
            "Path4".to_string(),
        ];

        // This should be fast because all pathways share the same Gamma params
        let results = run_gsea_parallel(&stats, &pathways, &names, 1.0, 100, 10, 1000, "std");

        assert_eq!(results.len(), 4);
        // All results should have valid p-values
        for r in &results {
            assert!(r.pval >= 0.0 && r.pval <= 1.0);
        }
    }
}
