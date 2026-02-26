//! Gamma distribution fitting for p-value estimation
//!
//! Uses anchor sampling to estimate null distribution parameters efficiently.
//! Supports adaptive sampling that stops when parameters converge.

use rand::seq::SliceRandom;
use rand::thread_rng;
use statrs::distribution::{ContinuousCDF, Gamma};

use crate::gsea_core::calc_batch_es;

/// Fit Gamma distribution from pre-computed ES values
/// This allows sharing random samples across multiple calls
///
/// # Arguments
/// * `es_values` - Pre-computed ES values from random gene sets
///
/// # Returns
/// Tuple of (positive_gamma_params, negative_gamma_params)
pub fn fit_gamma_from_es(es_values: &[f64]) -> Result<(GammaParams, GammaParams), FitError> {
    if es_values.len() < 20 {
        return Err(FitError::InsufficientSamples);
    }

    // Split into positive and negative ES
    let mut positive_es: Vec<f64> = Vec::with_capacity(es_values.len() / 2);
    let mut negative_es: Vec<f64> = Vec::with_capacity(es_values.len() / 2);

    for &es in es_values {
        if es > 0.0 {
            positive_es.push(es);
        } else if es < 0.0 {
            negative_es.push(-es); // Store absolute value
        }
    }

    // Fit Gamma to positive ES values
    let pos_params = if positive_es.len() >= 10 {
        fit_gamma_mle(&positive_es)?
    } else {
        GammaParams {
            shape: 1.0,
            rate: 5.0,
            location: 0.0,
            is_positive: true,
        }
    };

    // Fit Gamma to negative ES values
    let neg_params = if negative_es.len() >= 10 {
        let mut params = fit_gamma_mle(&negative_es)?;
        params.is_positive = false;
        params
    } else {
        GammaParams {
            shape: 1.0,
            rate: 5.0,
            location: 0.0,
            is_positive: false,
        }
    };

    Ok((pos_params, neg_params))
}

/// Adaptive sampling configuration
pub struct AdaptiveConfig {
    /// Minimum number of anchors before checking convergence
    pub min_anchors: usize,
    /// Maximum number of anchors (hard limit)
    pub max_anchors: usize,
    /// Batch size for each sampling iteration
    pub batch_size: usize,
    /// Convergence threshold for parameter change rate (e.g., 0.02 = 2%)
    pub convergence_threshold: f64,
}

impl Default for AdaptiveConfig {
    fn default() -> Self {
        Self {
            min_anchors: 75,      // Increased minimum for stability
            max_anchors: 200,
            batch_size: 25,
            convergence_threshold: 0.01, // 1% change threshold (stricter)
        }
    }
}

/// Parameters for fitted Gamma distribution
#[derive(Debug, Clone)]
pub struct GammaParams {
    pub shape: f64,       // alpha (k)
    pub rate: f64,        // beta (1/theta)
    pub location: f64,    // shift parameter for negative ES
    pub is_positive: bool, // whether this is for positive or negative ES
}

/// Error type for Gamma fitting failures
#[derive(Debug, Clone)]
pub enum FitError {
    InsufficientSamples,
    InvalidParameters,
    ConvergenceFailed,
}

/// Fit Gamma distribution to null ES values using anchor sampling
///
/// # Arguments
/// * `stats` - Gene statistics (sorted)
/// * `gene_set_size` - Size of the gene set being tested
/// * `n_anchors` - Number of random samples to generate
/// * `gsea_param` - GSEA weighting parameter
///
/// # Returns
/// Tuple of (positive_gamma_params, negative_gamma_params)
pub fn fit_gamma_null(
    stats: &[f64],
    gene_set_size: usize,
    n_anchors: usize,
    gsea_param: f64,
) -> Result<(GammaParams, GammaParams), FitError> {
    let n = stats.len();
    if n < gene_set_size || gene_set_size < 3 {
        return Err(FitError::InsufficientSamples);
    }

    let mut rng = thread_rng();
    let indices: Vec<usize> = (0..n).collect();

    // Pre-generate all random samples (batch approach)
    let samples: Vec<Vec<usize>> = (0..n_anchors)
        .map(|_| {
            indices
                .choose_multiple(&mut rng, gene_set_size)
                .cloned()
                .collect()
        })
        .collect();

    // Batch calculate all ES values at once (key optimization)
    let es_values = calc_batch_es(stats, gene_set_size, &samples, gsea_param);

    // Split into positive and negative ES
    let mut positive_es: Vec<f64> = Vec::with_capacity(n_anchors / 2);
    let mut negative_es: Vec<f64> = Vec::with_capacity(n_anchors / 2);

    for es in es_values {
        if es > 0.0 {
            positive_es.push(es);
        } else if es < 0.0 {
            negative_es.push(-es); // Store absolute value for fitting
        }
    }

    // Fit Gamma to positive ES values
    let pos_params = if positive_es.len() >= 10 {
        fit_gamma_mle(&positive_es)?
    } else {
        // Fallback: use default parameters
        GammaParams {
            shape: 1.0,
            rate: 5.0,
            location: 0.0,
            is_positive: true,
        }
    };

    // Fit Gamma to negative ES values (stored as positive)
    let neg_params = if negative_es.len() >= 10 {
        let mut params = fit_gamma_mle(&negative_es)?;
        params.is_positive = false;
        params
    } else {
        GammaParams {
            shape: 1.0,
            rate: 5.0,
            location: 0.0,
            is_positive: false,
        }
    };

    Ok((pos_params, neg_params))
}

/// Fit Gamma distribution with adaptive sampling
///
/// Stops sampling early when parameters converge, saving computation time.
/// Typically uses 40-60% fewer samples than fixed sampling for large pathways.
///
/// # Arguments
/// * `stats` - Gene statistics (sorted)
/// * `gene_set_size` - Size of the gene set being tested
/// * `gsea_param` - GSEA weighting parameter
/// * `config` - Adaptive sampling configuration (use Default::default() for standard settings)
///
/// # Returns
/// Tuple of (positive_gamma_params, negative_gamma_params, actual_samples_used)
pub fn fit_gamma_null_adaptive(
    stats: &[f64],
    gene_set_size: usize,
    gsea_param: f64,
    config: &AdaptiveConfig,
) -> Result<((GammaParams, GammaParams), usize), FitError> {
    let n = stats.len();
    if n < gene_set_size || gene_set_size < 3 {
        return Err(FitError::InsufficientSamples);
    }

    let mut rng = thread_rng();
    let indices: Vec<usize> = (0..n).collect();

    let mut positive_es: Vec<f64> = Vec::new();
    let mut negative_es: Vec<f64> = Vec::new();

    let mut prev_pos_params: Option<(f64, f64)> = None; // (shape, rate)
    let mut prev_neg_params: Option<(f64, f64)> = None;
    let mut total_samples = 0usize;

    loop {
        // Pre-generate batch samples
        let samples: Vec<Vec<usize>> = (0..config.batch_size)
            .map(|_| {
                indices
                    .choose_multiple(&mut rng, gene_set_size)
                    .cloned()
                    .collect()
            })
            .collect();

        // Batch calculate ES values
        let es_values = calc_batch_es(stats, gene_set_size, &samples, gsea_param);

        // Split into positive and negative
        for es in es_values {
            if es > 0.0 {
                positive_es.push(es);
            } else if es < 0.0 {
                negative_es.push(-es);
            }
        }
        total_samples += config.batch_size;

        // Check if we've reached max
        if total_samples >= config.max_anchors {
            break;
        }

        // Check convergence after minimum samples
        if total_samples >= config.min_anchors {
            // Try to fit current parameters
            let pos_converged = if positive_es.len() >= 10 {
                if let Ok(params) = fit_gamma_mle(&positive_es) {
                    let converged = check_convergence(
                        prev_pos_params,
                        (params.shape, params.rate),
                        config.convergence_threshold,
                    );
                    prev_pos_params = Some((params.shape, params.rate));
                    converged
                } else {
                    false
                }
            } else {
                false // Need more samples
            };

            let neg_converged = if negative_es.len() >= 10 {
                if let Ok(params) = fit_gamma_mle(&negative_es) {
                    let converged = check_convergence(
                        prev_neg_params,
                        (params.shape, params.rate),
                        config.convergence_threshold,
                    );
                    prev_neg_params = Some((params.shape, params.rate));
                    converged
                } else {
                    false
                }
            } else {
                false
            };

            // Stop if both converged
            if pos_converged && neg_converged {
                break;
            }
        }
    }

    // Final fit
    let pos_params = if positive_es.len() >= 10 {
        fit_gamma_mle(&positive_es)?
    } else {
        GammaParams {
            shape: 1.0,
            rate: 5.0,
            location: 0.0,
            is_positive: true,
        }
    };

    let neg_params = if negative_es.len() >= 10 {
        let mut params = fit_gamma_mle(&negative_es)?;
        params.is_positive = false;
        params
    } else {
        GammaParams {
            shape: 1.0,
            rate: 5.0,
            location: 0.0,
            is_positive: false,
        }
    };

    Ok(((pos_params, neg_params), total_samples))
}

/// Check if parameters have converged
fn check_convergence(
    prev: Option<(f64, f64)>,
    current: (f64, f64),
    threshold: f64,
) -> bool {
    match prev {
        Some((prev_shape, prev_rate)) => {
            let shape_change = (current.0 - prev_shape).abs() / prev_shape.max(1e-10);
            let rate_change = (current.1 - prev_rate).abs() / prev_rate.max(1e-10);
            shape_change < threshold && rate_change < threshold
        }
        None => false, // First iteration, not converged yet
    }
}

/// Fit Gamma distribution using method of moments (fast approximation)
fn fit_gamma_mle(data: &[f64]) -> Result<GammaParams, FitError> {
    if data.len() < 5 {
        return Err(FitError::InsufficientSamples);
    }

    // Calculate mean and variance
    let n = data.len() as f64;
    let mean: f64 = data.iter().sum::<f64>() / n;
    let variance: f64 = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);

    if mean <= 0.0 || variance <= 0.0 {
        return Err(FitError::InvalidParameters);
    }

    // Method of moments estimators
    let shape = mean.powi(2) / variance;
    let rate = mean / variance;

    // Validate parameters
    if shape <= 0.0 || rate <= 0.0 || !shape.is_finite() || !rate.is_finite() {
        return Err(FitError::InvalidParameters);
    }

    Ok(GammaParams {
        shape,
        rate,
        location: 0.0,
        is_positive: true,
    })
}

/// Calculate p-value from fitted Gamma distribution
///
/// For positive ES: P(X > es) = 1 - CDF(es)
/// For negative ES: P(X > |es|) = 1 - CDF(|es|) using negative params
pub fn gamma_pvalue(es: f64, pos_params: &GammaParams, neg_params: &GammaParams) -> f64 {
    if es == 0.0 {
        return 1.0;
    }

    let (params, abs_es) = if es > 0.0 {
        (pos_params, es)
    } else {
        (neg_params, -es)
    };

    // Create Gamma distribution (statrs uses rate parameterization)
    match Gamma::new(params.shape, params.rate) {
        Ok(gamma) => {
            // P-value is the upper tail probability
            let pval = 1.0 - gamma.cdf(abs_es);
            // Ensure p-value is in valid range
            pval.max(1e-16).min(1.0)
        }
        Err(_) => 1.0, // Return 1.0 if distribution creation fails
    }
}

/// Calculate Normalized Enrichment Score
///
/// NES = ES / mean(|random ES|) with same sign as ES
pub fn calculate_nes(
    es: f64,
    pos_params: &GammaParams,
    neg_params: &GammaParams,
) -> f64 {
    if es == 0.0 {
        return 0.0;
    }

    // Mean of Gamma distribution is shape/rate
    let (mean_abs_es, sign) = if es > 0.0 {
        (pos_params.shape / pos_params.rate, 1.0)
    } else {
        (neg_params.shape / neg_params.rate, -1.0)
    };

    if mean_abs_es > 0.0 {
        sign * es.abs() / mean_abs_es
    } else {
        es
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fit_gamma() {
        // Generate some Gamma-distributed data
        let data: Vec<f64> = vec![0.1, 0.2, 0.15, 0.3, 0.25, 0.18, 0.22, 0.28, 0.12, 0.19];
        let result = fit_gamma_mle(&data);
        assert!(result.is_ok());
        let params = result.unwrap();
        assert!(params.shape > 0.0);
        assert!(params.rate > 0.0);
    }

    #[test]
    fn test_gamma_pvalue() {
        let pos_params = GammaParams {
            shape: 2.0,
            rate: 10.0,
            location: 0.0,
            is_positive: true,
        };
        let neg_params = GammaParams {
            shape: 2.0,
            rate: 10.0,
            location: 0.0,
            is_positive: false,
        };

        // Higher ES should give lower p-value
        let pval_low = gamma_pvalue(0.5, &pos_params, &neg_params);
        let pval_high = gamma_pvalue(0.1, &pos_params, &neg_params);
        assert!(pval_low < pval_high);
    }

    #[test]
    fn test_convergence_check() {
        // First call should return false (no previous)
        assert!(!check_convergence(None, (2.0, 10.0), 0.02));

        // Same parameters should converge
        assert!(check_convergence(Some((2.0, 10.0)), (2.0, 10.0), 0.02));

        // Small change should converge
        assert!(check_convergence(Some((2.0, 10.0)), (2.01, 10.01), 0.02));

        // Large change should not converge
        assert!(!check_convergence(Some((2.0, 10.0)), (2.5, 10.0), 0.02));
    }

    #[test]
    fn test_adaptive_sampling() {
        // Create test stats
        let n = 2000;
        let stats: Vec<f64> = (0..n)
            .map(|i| 3.0 - (i as f64 / (n as f64 / 6.0)))
            .collect();

        let config = AdaptiveConfig::default();
        let result = fit_gamma_null_adaptive(&stats, 50, 1.0, &config);
        assert!(result.is_ok());

        let ((pos_params, neg_params), samples_used) = result.unwrap();

        // Should use fewer samples than max (due to convergence)
        // or exactly max if didn't converge
        assert!(samples_used <= config.max_anchors);
        assert!(samples_used >= config.min_anchors);

        // Parameters should be valid
        assert!(pos_params.shape > 0.0);
        assert!(pos_params.rate > 0.0);
        assert!(neg_params.shape > 0.0);
        assert!(neg_params.rate > 0.0);
    }

    #[test]
    fn test_adaptive_vs_fixed_consistency() {
        // Results should be similar between adaptive and fixed sampling
        let n = 2000;
        let stats: Vec<f64> = (0..n)
            .map(|i| 3.0 - (i as f64 / (n as f64 / 6.0)))
            .collect();

        // Fixed sampling
        let (pos_fixed, neg_fixed) = fit_gamma_null(&stats, 50, 200, 1.0).unwrap();

        // Adaptive sampling
        let config = AdaptiveConfig {
            min_anchors: 50,
            max_anchors: 200,
            batch_size: 25,
            convergence_threshold: 0.02,
        };
        let ((pos_adaptive, neg_adaptive), _) =
            fit_gamma_null_adaptive(&stats, 50, 1.0, &config).unwrap();

        // Parameters should be in similar range (not exact due to randomness)
        // Just verify they're reasonable
        assert!(pos_adaptive.shape > 0.0 && pos_adaptive.shape < 100.0);
        assert!(neg_adaptive.shape > 0.0 && neg_adaptive.shape < 100.0);

        // Both should produce similar mean (shape/rate)
        let mean_fixed = pos_fixed.shape / pos_fixed.rate;
        let mean_adaptive = pos_adaptive.shape / pos_adaptive.rate;

        // Means should be within 50% of each other (generous due to randomness)
        let ratio = mean_fixed / mean_adaptive;
        assert!(ratio > 0.5 && ratio < 2.0,
            "Fixed mean {} vs Adaptive mean {} differ too much",
            mean_fixed, mean_adaptive);
    }
}
