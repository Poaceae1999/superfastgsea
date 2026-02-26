//! Empirical CDF + GPD (Generalized Pareto Distribution) for p-value estimation
//!
//! This method uses:
//! - Empirical CDF for the center of the distribution (80%)
//! - GPD fitting for tails (10% each side) based on Extreme Value Theory
//!
//! Advantages over pure Gamma:
//! - No distributional assumption for the center
//! - EVT provides theoretical guarantees for tail estimation
//! - More accurate p-values for extreme enrichment scores

use rand::seq::SliceRandom;
use rand::thread_rng;

use crate::gsea_core::calc_enrichment_score_sparse;

/// Parameters for fitted GPD (Generalized Pareto Distribution)
/// F(x) = 1 - (1 + ξx/σ)^(-1/ξ) for ξ ≠ 0
/// F(x) = 1 - exp(-x/σ) for ξ = 0
#[derive(Debug, Clone)]
pub struct GpdParams {
    pub shape: f64,     // ξ (xi) - shape parameter
    pub scale: f64,     // σ (sigma) - scale parameter
    pub threshold: f64, // u - threshold for exceedances
    pub n_exceedances: usize, // number of samples above threshold
}

/// Combined null distribution parameters
#[derive(Debug, Clone)]
pub struct EmpiricalGpdParams {
    pub positive_samples: Vec<f64>,  // sorted positive ES samples
    pub negative_samples: Vec<f64>,  // sorted absolute negative ES samples
    pub pos_gpd: Option<GpdParams>,  // GPD for positive tail
    pub neg_gpd: Option<GpdParams>,  // GPD for negative tail
    pub pos_threshold_quantile: f64, // quantile used for positive tail threshold
    pub neg_threshold_quantile: f64, // quantile used for negative tail threshold
}

/// Error type for GPD fitting
#[derive(Debug, Clone)]
pub enum GpdFitError {
    InsufficientSamples,
    InvalidParameters,
    NumericalError,
}

/// Minimum samples needed for reliable GPD estimation
const MIN_TAIL_SAMPLES: usize = 30;

/// Default tail proportion (10% each side)
const DEFAULT_TAIL_PROPORTION: f64 = 0.10;

/// Fit GPD to exceedances using Probability Weighted Moments (PWM) method
///
/// PWM estimators are more robust than MLE for small samples
fn fit_gpd_pwm(exceedances: &[f64]) -> Result<GpdParams, GpdFitError> {
    let n = exceedances.len();
    if n < MIN_TAIL_SAMPLES {
        return Err(GpdFitError::InsufficientSamples);
    }

    // Sort exceedances
    let mut sorted: Vec<f64> = exceedances.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Calculate probability weighted moments
    // β_0 = mean
    // β_1 = (1/n) Σ [(i-1)/(n-1)] * x_(i)
    let beta_0: f64 = sorted.iter().sum::<f64>() / n as f64;

    let beta_1: f64 = sorted
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let p = i as f64 / (n - 1) as f64;
            p * x
        })
        .sum::<f64>() / n as f64;

    // PWM estimators for GPD
    // ξ = 2 - β_0 / (β_0 - 2β_1)
    // σ = 2β_0β_1 / (β_0 - 2β_1)

    let denom = beta_0 - 2.0 * beta_1;

    if denom.abs() < 1e-10 || denom <= 0.0 {
        // Fallback to exponential (ξ = 0)
        return Ok(GpdParams {
            shape: 0.0,
            scale: beta_0,
            threshold: 0.0,
            n_exceedances: n,
        });
    }

    let shape = 2.0 - beta_0 / denom;
    let scale = 2.0 * beta_0 * beta_1 / denom;

    // Validate parameters
    if !scale.is_finite() || scale <= 0.0 {
        return Err(GpdFitError::InvalidParameters);
    }

    // Bound shape parameter to reasonable range [-0.5, 0.5]
    // Most ES distributions have light to moderate tails
    let bounded_shape = shape.max(-0.5).min(0.5);

    Ok(GpdParams {
        shape: bounded_shape,
        scale,
        threshold: 0.0, // Will be set by caller
        n_exceedances: n,
    })
}

/// Calculate GPD survival probability P(X > x)
fn gpd_survival(x: f64, params: &GpdParams) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }

    let xi = params.shape;
    let sigma = params.scale;

    if xi.abs() < 1e-10 {
        // Exponential case (ξ → 0)
        (-x / sigma).exp()
    } else {
        let arg = 1.0 + xi * x / sigma;
        if arg <= 0.0 {
            // Beyond support of distribution
            0.0
        } else {
            arg.powf(-1.0 / xi)
        }
    }
}

/// Fit empirical distribution + GPD for null ES values
///
/// # Arguments
/// * `stats` - Gene statistics (sorted)
/// * `gene_set_size` - Size of the gene set being tested
/// * `n_anchors` - Number of random samples (recommended ≥ 500)
/// * `gsea_param` - GSEA weighting parameter
/// * `tail_proportion` - Proportion of samples for tail fitting (default 0.10)
///
/// # Returns
/// Combined empirical + GPD parameters
pub fn fit_empirical_gpd_null(
    stats: &[f64],
    gene_set_size: usize,
    n_anchors: usize,
    gsea_param: f64,
    tail_proportion: Option<f64>,
) -> Result<EmpiricalGpdParams, GpdFitError> {
    let n = stats.len();
    if n < gene_set_size || gene_set_size < 3 {
        return Err(GpdFitError::InsufficientSamples);
    }

    let tail_prop = tail_proportion.unwrap_or(DEFAULT_TAIL_PROPORTION);

    let mut rng = thread_rng();
    let indices: Vec<usize> = (0..n).collect();

    // Generate random ES values
    let mut positive_es: Vec<f64> = Vec::new();
    let mut negative_es: Vec<f64> = Vec::new();

    for _ in 0..n_anchors {
        let sample: Vec<usize> = indices
            .choose_multiple(&mut rng, gene_set_size)
            .cloned()
            .collect();

        let (es, _) = calc_enrichment_score_sparse(stats, &sample, gsea_param);

        if es > 0.0 {
            positive_es.push(es);
        } else if es < 0.0 {
            negative_es.push(-es); // Store absolute value
        }
    }

    // Sort samples for empirical CDF
    positive_es.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    negative_es.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Fit GPD for positive tail
    let pos_gpd = if positive_es.len() >= MIN_TAIL_SAMPLES {
        let threshold_idx = ((1.0 - tail_prop) * positive_es.len() as f64) as usize;
        let threshold_idx = threshold_idx.min(positive_es.len() - MIN_TAIL_SAMPLES);
        let threshold = positive_es[threshold_idx];

        let exceedances: Vec<f64> = positive_es[threshold_idx..]
            .iter()
            .map(|&x| x - threshold)
            .collect();

        match fit_gpd_pwm(&exceedances) {
            Ok(mut params) => {
                params.threshold = threshold;
                Some(params)
            }
            Err(_) => None,
        }
    } else {
        None
    };

    // Fit GPD for negative tail
    let neg_gpd = if negative_es.len() >= MIN_TAIL_SAMPLES {
        let threshold_idx = ((1.0 - tail_prop) * negative_es.len() as f64) as usize;
        let threshold_idx = threshold_idx.min(negative_es.len() - MIN_TAIL_SAMPLES);
        let threshold = negative_es[threshold_idx];

        let exceedances: Vec<f64> = negative_es[threshold_idx..]
            .iter()
            .map(|&x| x - threshold)
            .collect();

        match fit_gpd_pwm(&exceedances) {
            Ok(mut params) => {
                params.threshold = threshold;
                Some(params)
            }
            Err(_) => None,
        }
    } else {
        None
    };

    Ok(EmpiricalGpdParams {
        positive_samples: positive_es,
        negative_samples: negative_es,
        pos_gpd,
        neg_gpd,
        pos_threshold_quantile: 1.0 - tail_prop,
        neg_threshold_quantile: 1.0 - tail_prop,
    })
}

/// Calculate p-value using empirical CDF + GPD
///
/// For values in the center: use empirical CDF directly
/// For values in the tail: use GPD extrapolation
pub fn empirical_gpd_pvalue(es: f64, params: &EmpiricalGpdParams) -> f64 {
    if es == 0.0 {
        return 1.0;
    }

    let (samples, gpd_params, abs_es) = if es > 0.0 {
        (&params.positive_samples, &params.pos_gpd, es)
    } else {
        (&params.negative_samples, &params.neg_gpd, -es)
    };

    if samples.is_empty() {
        return 1.0;
    }

    let n = samples.len();

    // Find position in sorted samples using binary search
    let rank = match samples.binary_search_by(|x| {
        x.partial_cmp(&abs_es).unwrap_or(std::cmp::Ordering::Equal)
    }) {
        Ok(pos) => pos,
        Err(pos) => pos,
    };

    // Empirical survival probability
    let empirical_p = (n - rank) as f64 / n as f64;

    // If we have GPD params and value exceeds threshold, use GPD for tail
    if let Some(gpd) = gpd_params {
        if abs_es > gpd.threshold {
            // P(X > x) = P(X > u) * P(X > x | X > u)
            //          = tail_proportion * GPD_survival(x - u)
            let tail_prob = gpd.n_exceedances as f64 / n as f64;
            let exceedance = abs_es - gpd.threshold;
            let gpd_surv = gpd_survival(exceedance, gpd);
            let pval = tail_prob * gpd_surv;
            return pval.max(1e-16).min(1.0);
        }
    }

    // Use empirical CDF for center
    empirical_p.max(1e-16).min(1.0)
}

/// Calculate NES using empirical distribution
///
/// NES = ES / mean(|null ES|) with same sign
pub fn calculate_nes_empirical(es: f64, params: &EmpiricalGpdParams) -> f64 {
    if es == 0.0 {
        return 0.0;
    }

    let (samples, sign) = if es > 0.0 {
        (&params.positive_samples, 1.0)
    } else {
        (&params.negative_samples, -1.0)
    };

    if samples.is_empty() {
        return es;
    }

    let mean_abs_es: f64 = samples.iter().sum::<f64>() / samples.len() as f64;

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
    fn test_gpd_pwm_fit() {
        // Generate exponential-like exceedances (GPD with ξ ≈ 0)
        let exceedances: Vec<f64> = (0..100)
            .map(|i| 0.01 * (i as f64 + 1.0))
            .collect();

        let result = fit_gpd_pwm(&exceedances);
        assert!(result.is_ok());

        let params = result.unwrap();
        assert!(params.scale > 0.0);
        // Shape should be close to 0 for exponential data
        assert!(params.shape.abs() < 0.5);
    }

    #[test]
    fn test_gpd_survival() {
        // Exponential case (ξ = 0)
        let params = GpdParams {
            shape: 0.0,
            scale: 0.2,
            threshold: 0.0,
            n_exceedances: 100,
        };

        // P(X > 0) should be 1
        assert!((gpd_survival(0.0, &params) - 1.0).abs() < 1e-10);

        // Survival should decrease with x
        let p1 = gpd_survival(0.1, &params);
        let p2 = gpd_survival(0.2, &params);
        assert!(p1 > p2);
        assert!(p1 > 0.0 && p1 < 1.0);
    }

    #[test]
    fn test_empirical_gpd_pvalue_ordering() {
        // Create mock params
        let positive_samples: Vec<f64> = (1..=100).map(|i| i as f64 * 0.01).collect();
        let negative_samples: Vec<f64> = (1..=100).map(|i| i as f64 * 0.01).collect();

        let params = EmpiricalGpdParams {
            positive_samples,
            negative_samples,
            pos_gpd: None,
            neg_gpd: None,
            pos_threshold_quantile: 0.9,
            neg_threshold_quantile: 0.9,
        };

        // Higher ES should give lower p-value
        let p1 = empirical_gpd_pvalue(0.5, &params);
        let p2 = empirical_gpd_pvalue(0.8, &params);
        assert!(p2 < p1, "Higher ES should have lower p-value");
    }

    #[test]
    fn test_nes_calculation() {
        let positive_samples: Vec<f64> = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        let negative_samples: Vec<f64> = vec![0.1, 0.2, 0.3];

        let params = EmpiricalGpdParams {
            positive_samples,
            negative_samples,
            pos_gpd: None,
            neg_gpd: None,
            pos_threshold_quantile: 0.9,
            neg_threshold_quantile: 0.9,
        };

        // Mean of positive samples is 0.3
        let nes = calculate_nes_empirical(0.6, &params);
        assert!((nes - 2.0).abs() < 1e-10, "NES should be 0.6/0.3 = 2.0");

        // Negative ES
        let nes_neg = calculate_nes_empirical(-0.4, &params);
        assert!(nes_neg < 0.0, "NES should be negative for negative ES");
    }
}
