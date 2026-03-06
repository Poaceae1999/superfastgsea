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

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// Cached parameters for either method
enum CachedParams {
    Gamma(Option<(GammaParams, GammaParams)>),
    EmpiricalGpd(Option<EmpiricalGpdParams>),
}

/// Get the representative size for a pathway size (bin center)
#[inline]
fn get_size_bin(size: usize) -> usize {
    let bin = size / SIZE_BIN_WIDTH;
    bin * SIZE_BIN_WIDTH + SIZE_BIN_WIDTH / 2
}

/// Filter pathways by size bounds, returning (original_index, pathway, name) triples.
fn filter_pathways_by_size<'a>(
    pathways: &'a [Vec<usize>],
    pathway_names: &'a [String],
    min_size: usize,
    max_size: usize,
) -> Vec<(usize, &'a Vec<usize>, &'a String)> {
    pathways
        .iter()
        .zip(pathway_names.iter())
        .enumerate()
        .filter(|(_, (pathway, _))| {
            let size = pathway.len();
            size >= min_size && size <= max_size
        })
        .map(|(i, (pathway, name))| (i, pathway, name))
        .collect()
}

/// Collect unique size bins from filtered pathways (sorted, deduplicated,
/// excluding bins below MIN_SIZE_FOR_GAMMA).
fn collect_unique_size_bins(filtered: &[(usize, &Vec<usize>, &String)]) -> Vec<usize> {
    let mut bins: Vec<usize> = filtered
        .iter()
        .map(|(_, pathway, _)| get_size_bin(pathway.len()))
        .filter(|&bin| bin >= MIN_SIZE_FOR_GAMMA)
        .collect();
    bins.sort_unstable();
    bins.dedup();
    bins
}

/// Apply score type filter to an enrichment score.
#[inline]
fn apply_score_type(es: f64, score_type: &str) -> f64 {
    match score_type {
        "pos" if es < 0.0 => 0.0,
        "neg" if es > 0.0 => 0.0,
        _ => es,
    }
}

/// Build a Gamma parameter cache for the given size bins using fixed-anchor fitting.
fn build_gamma_cache(
    stats: &[f64],
    unique_bins: &[usize],
    n_anchors: usize,
    gsea_param: f64,
) -> HashMap<usize, Option<(GammaParams, GammaParams)>> {
    unique_bins
        .par_iter()
        .map(|&bin_size| {
            let params = fit_gamma_null(stats, bin_size, n_anchors, gsea_param).ok();
            (bin_size, params)
        })
        .collect()
}

/// Compute p-value and NES for a single pathway using a method-aware cache.
/// Falls back to permutation when the gene set is too small or fitting failed.
fn compute_pval_nes(
    stats: &[f64],
    pathway: &[usize],
    filtered_es: f64,
    gsea_param: f64,
    cache: &HashMap<usize, CachedParams>,
) -> (f64, f64) {
    let size = pathway.len();
    if size < MIN_SIZE_FOR_GAMMA {
        return permutation_fallback(stats, pathway, filtered_es, gsea_param);
    }
    let lookup = get_size_bin(size);
    match cache.get(&lookup) {
        Some(CachedParams::Gamma(Some((pos, neg)))) => {
            (gamma_pvalue(filtered_es, pos, neg), calculate_nes(filtered_es, pos, neg))
        }
        Some(CachedParams::EmpiricalGpd(Some(params))) => {
            (empirical_gpd_pvalue(filtered_es, params), calculate_nes_empirical(filtered_es, params))
        }
        _ => permutation_fallback(stats, pathway, filtered_es, gsea_param),
    }
}

/// Compute p-value and NES for a single pathway using a Gamma-only cache.
/// Falls back to permutation when the gene set is too small or fitting failed.
fn compute_pval_nes_gamma(
    stats: &[f64],
    pathway: &[usize],
    filtered_es: f64,
    gsea_param: f64,
    gamma_cache: &HashMap<usize, Option<(GammaParams, GammaParams)>>,
) -> (f64, f64) {
    let size = pathway.len();
    if size < MIN_SIZE_FOR_GAMMA {
        return permutation_fallback(stats, pathway, filtered_es, gsea_param);
    }
    let bin_size = get_size_bin(size);
    match gamma_cache.get(&bin_size) {
        Some(Some((pos, neg))) => {
            (gamma_pvalue(filtered_es, pos, neg), calculate_nes(filtered_es, pos, neg))
        }
        _ => permutation_fallback(stats, pathway, filtered_es, gsea_param),
    }
}

/// Permutation-based p-value/NES fallback for small or failed-fit gene sets.
#[inline]
fn permutation_fallback(stats: &[f64], pathway: &[usize], filtered_es: f64, gsea_param: f64) -> (f64, f64) {
    let (p, mean_pos, mean_neg) = permutation_pvalue(stats, pathway, gsea_param, DEFAULT_N_PERM);
    let n = calculate_nes_permutation(filtered_es, mean_pos, mean_neg);
    (p, n)
}

/// Analyze a single pathway using a Gamma-only cache.
fn analyze_pathway_gamma(
    stats: &[f64],
    pathway: &[usize],
    name: &str,
    gsea_param: f64,
    score_type: &str,
    gamma_cache: &HashMap<usize, Option<(GammaParams, GammaParams)>>,
) -> GseaResult {
    let (es, leading_edge) = calc_enrichment_score_sparse(stats, pathway, gsea_param);
    let filtered_es = apply_score_type(es, score_type);
    let (pval, nes) = compute_pval_nes_gamma(stats, pathway, filtered_es, gsea_param, gamma_cache);

    GseaResult {
        pathway: name.to_string(),
        es: filtered_es,
        nes,
        pval,
        size: pathway.len(),
        leading_edge,
    }
}

/// Analyze a single pathway using a method-aware cache.
fn analyze_pathway_with_cache(
    stats: &[f64],
    pathway: &[usize],
    name: &str,
    gsea_param: f64,
    score_type: &str,
    cache: &HashMap<usize, CachedParams>,
) -> GseaResult {
    let (es, leading_edge) = calc_enrichment_score_sparse(stats, pathway, gsea_param);
    let filtered_es = apply_score_type(es, score_type);
    let (pval, nes) = compute_pval_nes(stats, pathway, filtered_es, gsea_param, cache);

    GseaResult {
        pathway: name.to_string(),
        es: filtered_es,
        nes,
        pval,
        size: pathway.len(),
        leading_edge,
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Run GSEA analysis on multiple pathways in parallel with Gamma parameter caching
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
    let effective_anchors = n_anchors.max(DEFAULT_N_ANCHORS);
    let filtered = filter_pathways_by_size(pathways, pathway_names, min_size, max_size);
    let unique_bins = collect_unique_size_bins(&filtered);
    let gamma_cache = build_gamma_cache(stats, &unique_bins, effective_anchors, gsea_param);

    filtered
        .par_iter()
        .map(|(_, pathway, name)| {
            analyze_pathway_gamma(stats, pathway, name, gsea_param, score_type, &gamma_cache)
        })
        .collect()
}

/// Run GSEA analysis with method selection (Gamma or Empirical+GPD)
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
    let effective_anchors = match method {
        PvalueMethod::Gamma => n_anchors.max(DEFAULT_N_ANCHORS),
        PvalueMethod::EmpiricalGpd => n_anchors.max(DEFAULT_N_ANCHORS_GPD),
    };

    let filtered = filter_pathways_by_size(pathways, pathway_names, min_size, max_size);
    let unique_bins = collect_unique_size_bins(&filtered);

    let params_cache: HashMap<usize, CachedParams> = match method {
        PvalueMethod::Gamma => {
            unique_bins
                .par_iter()
                .map(|&bin_size| {
                    let params = fit_gamma_null(stats, bin_size, effective_anchors, gsea_param).ok();
                    (bin_size, CachedParams::Gamma(params))
                })
                .collect()
        }
        PvalueMethod::EmpiricalGpd => {
            unique_bins
                .par_iter()
                .map(|&bin_size| {
                    let params = fit_empirical_gpd_null(
                        stats, bin_size, effective_anchors, gsea_param, None,
                    )
                    .ok();
                    (bin_size, CachedParams::EmpiricalGpd(params))
                })
                .collect()
        }
    };

    filtered
        .par_iter()
        .map(|(_, pathway, name)| {
            analyze_pathway_with_cache(stats, pathway, name, gsea_param, score_type, &params_cache)
        })
        .collect()
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
    let filtered = filter_pathways_by_size(pathways, pathway_names, min_size, max_size);
    if filtered.is_empty() {
        return (Vec::new(), 0, 0);
    }

    let unique_bins = collect_unique_size_bins(&filtered);

    // === PASS 1: Coarse estimation ===
    let coarse_cache = build_gamma_cache(stats, &unique_bins, config.coarse_anchors, gsea_param);

    let coarse_results: Vec<(usize, GseaResult, bool)> = filtered
        .par_iter()
        .map(|(idx, pathway, name)| {
            let result = analyze_pathway_gamma(
                stats, pathway, name, gsea_param, score_type, &coarse_cache,
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
            if let Some(&(_, pathway, name)) = filtered.iter().find(|(i, _, _)| *i == idx) {
                candidates.push((idx, pathway, name));
            }
        } else {
            final_results.push((idx, result));
        }
    }

    let coarse_count = filtered.len() - candidates.len();
    let fine_count = candidates.len();

    // === PASS 2: Fine estimation for candidates only ===
    if !candidates.is_empty() {
        let candidate_bins = collect_unique_size_bins(&candidates);
        let fine_cache = build_gamma_cache(stats, &candidate_bins, config.fine_anchors, gsea_param);

        let fine_results: Vec<(usize, GseaResult)> = candidates
            .par_iter()
            .map(|(idx, pathway, name)| {
                let result = analyze_pathway_gamma(
                    stats, pathway, name, gsea_param, score_type, &fine_cache,
                );
                (*idx, result)
            })
            .collect();

        final_results.extend(fine_results);
    }

    // Sort by original index to maintain order
    final_results.sort_by_key(|(idx, _)| *idx);
    let results = final_results.into_iter().map(|(_, r)| r).collect();

    (results, coarse_count, fine_count)
}

/// Run GSEA with adaptive sampling (Gamma method only)
///
/// Uses adaptive sampling that stops when parameters converge,
/// typically using 40-60% fewer samples than fixed sampling.
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
    let filtered = filter_pathways_by_size(pathways, pathway_names, min_size, max_size);
    let unique_bins = collect_unique_size_bins(&filtered);

    // Pre-compute Gamma parameters with adaptive sampling, tracking total samples
    let results_with_samples: Vec<(usize, Option<(GammaParams, GammaParams)>, usize)> = unique_bins
        .par_iter()
        .map(|&bin_size| {
            match fit_gamma_null_adaptive(stats, bin_size, gsea_param, &config) {
                Ok((params, samples)) => (bin_size, Some(params), samples),
                Err(_) => (bin_size, None, 0),
            }
        })
        .collect();

    let mut gamma_cache: HashMap<usize, Option<(GammaParams, GammaParams)>> = HashMap::new();
    let mut total_samples = 0usize;
    for (bin_size, params, samples) in results_with_samples {
        gamma_cache.insert(bin_size, params);
        total_samples += samples;
    }

    let results = filtered
        .par_iter()
        .map(|(_, pathway, name)| {
            analyze_pathway_gamma(stats, pathway, name, gsea_param, score_type, &gamma_cache)
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
    let filtered = filter_pathways_by_size(pathways, pathway_names, min_size, max_size);
    if filtered.is_empty() {
        return Vec::new();
    }

    let unique_bins = collect_unique_size_bins(&filtered);

    // === BATCH OPTIMIZATION: Generate all random samples and compute ES at once ===
    let mut rng = thread_rng();
    let indices: Vec<usize> = (0..n).collect();

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

    // Batch compute ES for all samples, then fit Gamma
    let gamma_cache: HashMap<usize, Option<(GammaParams, GammaParams)>> = unique_bins
        .par_iter()
        .map(|&bin_size| {
            if let Some(samples) = samples_per_bin.get(&bin_size) {
                let es_values = calc_batch_es(stats, bin_size, samples, gsea_param);
                let params = fit_gamma_from_es(&es_values).ok();
                (bin_size, params)
            } else {
                (bin_size, None)
            }
        })
        .collect();

    filtered
        .par_iter()
        .map(|(_, pathway, name)| {
            analyze_pathway_gamma(stats, pathway, name, gsea_param, score_type, &gamma_cache)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parallel_gsea() {
        let n = 2000;
        let stats: Vec<f64> = (0..n)
            .map(|i| 3.0 - (i as f64 / (n as f64 / 6.0)))
            .collect();

        let pathways: Vec<Vec<usize>> = vec![
            (0..50).collect(),
            (975..1025).collect(),
            (1950..2000).collect(),
        ];

        let names: Vec<String> = vec![
            "TopPathway".to_string(),
            "MiddlePathway".to_string(),
            "BottomPathway".to_string(),
        ];

        let results = run_gsea_parallel(
            &stats, &pathways, &names, 1.0, 200, 10, 1000, "std",
        );

        assert_eq!(results.len(), 3);

        let top = results.iter().find(|r| r.pathway == "TopPathway").unwrap();
        assert!(top.es > 0.0, "Top pathway ES should be positive: {}", top.es);
        assert!(top.pval <= 1.0, "Top pathway p-value should be valid: {}", top.pval);

        let bottom = results.iter().find(|r| r.pathway == "BottomPathway").unwrap();
        assert!(bottom.es < 0.0, "Bottom pathway ES should be negative: {}", bottom.es);

        let middle = results.iter().find(|r| r.pathway == "MiddlePathway").unwrap();
        assert!(
            middle.es.abs() < top.es.abs(),
            "Middle pathway ES {} should be smaller than top {}",
            middle.es.abs(),
            top.es.abs()
        );
    }

    #[test]
    fn test_gamma_caching() {
        let n = 1000;
        let stats: Vec<f64> = (0..n).map(|i| 2.0 - (i as f64 / 250.0)).collect();

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

        let results = run_gsea_parallel(&stats, &pathways, &names, 1.0, 100, 10, 1000, "std");

        assert_eq!(results.len(), 4);
        for r in &results {
            assert!(r.pval >= 0.0 && r.pval <= 1.0);
        }
    }
}
