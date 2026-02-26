use extendr_api::prelude::*;

mod gsea_core;
mod gamma_fit;
mod permutation;
mod parallel;
mod empirical_gpd;

pub use gsea_core::*;
pub use gamma_fit::*;
pub use permutation::*;
pub use parallel::*;
pub use empirical_gpd::*;

/// Test function to verify Rust-R bridge works
/// @export
#[extendr]
fn rust_hello() -> String {
    "Hello from Rust! superfastGSEA is ready.".to_string()
}

/// Calculate enrichment score for a single gene set
/// @param stats Named numeric vector of gene statistics (sorted by stat value)
/// @param gene_set_indices Integer vector of indices (1-based) for genes in the set
/// @param gsea_param Weighting parameter for gene statistics (default 1.0)
/// @return List containing ES, NES estimate, and leading edge indices
/// @export
#[extendr]
fn calc_es_rust(
    stats: &[f64],
    gene_set_indices: &[i32],
    gsea_param: f64,
) -> List {
    // Convert 1-based R indices to 0-based Rust indices
    let indices: Vec<usize> = gene_set_indices
        .iter()
        .map(|&i| (i - 1) as usize)
        .filter(|&i| i < stats.len())
        .collect();

    let (es, leading_edge) = gsea_core::calc_enrichment_score_sparse(stats, &indices, gsea_param);

    // Convert back to 1-based indices for R
    let leading_edge_r: Vec<i32> = leading_edge.iter().map(|&i| (i + 1) as i32).collect();

    list!(
        ES = es,
        leadingEdge = leading_edge_r
    )
}

/// Run GSEA analysis on multiple pathways in parallel
/// @param stats Named numeric vector of gene statistics
/// @param pathway_indices List of integer vectors, each containing gene indices for a pathway
/// @param pathway_names Character vector of pathway names
/// @param gsea_param Weighting parameter (default 1.0)
/// @param n_anchors Number of anchor samples for Gamma fitting (default 200)
/// @param min_size Minimum pathway size (default 1)
/// @param max_size Maximum pathway size (default Inf)
/// @param score_type Score type: "std", "pos", or "neg" (default "std")
/// @return Data frame with GSEA results
/// @export
#[extendr]
fn run_gsea_rust(
    stats: &[f64],
    pathway_indices: List,
    pathway_names: Vec<String>,
    gsea_param: f64,
    n_anchors: i32,
    min_size: i32,
    max_size: i32,
    score_type: &str,
) -> Robj {
    // Convert pathway indices from R (1-based) to Rust (0-based)
    let pathways: Vec<Vec<usize>> = pathway_indices
        .iter()
        .map(|(_, robj)| {
            if let Some(indices) = robj.as_integer_slice() {
                indices
                    .iter()
                    .map(|&i| (i - 1) as usize)
                    .filter(|&i| i < stats.len())
                    .collect()
            } else {
                Vec::new()
            }
        })
        .collect();

    // Filter pathways by size
    let min_sz = min_size as usize;
    let max_sz = if max_size < 0 { usize::MAX } else { max_size as usize };

    let results = parallel::run_gsea_parallel(
        stats,
        &pathways,
        &pathway_names,
        gsea_param,
        n_anchors as usize,
        min_sz,
        max_sz,
        score_type,
    );

    // Build output vectors
    let n = results.len();
    let mut out_pathway: Vec<String> = Vec::with_capacity(n);
    let mut out_pval: Vec<f64> = Vec::with_capacity(n);
    let mut out_es: Vec<f64> = Vec::with_capacity(n);
    let mut out_nes: Vec<f64> = Vec::with_capacity(n);
    let mut out_size: Vec<i32> = Vec::with_capacity(n);
    let mut out_leading_edge: Vec<Robj> = Vec::with_capacity(n);

    for result in results {
        out_pathway.push(result.pathway);
        out_pval.push(result.pval);
        out_es.push(result.es);
        out_nes.push(result.nes);
        out_size.push(result.size as i32);
        // Convert leading edge back to 1-based
        let le: Vec<i32> = result.leading_edge.iter().map(|&i| (i + 1) as i32).collect();
        out_leading_edge.push(le.into_robj());
    }

    // Return as a list - R will handle data frame conversion
    list!(
        pathway = out_pathway,
        pval = out_pval,
        ES = out_es,
        NES = out_nes,
        size = out_size,
        leadingEdge = List::from_values(out_leading_edge)
    ).into()
}

/// Run GSEA analysis with method selection
/// @param stats Named numeric vector of gene statistics
/// @param pathway_indices List of integer vectors, each containing gene indices for a pathway
/// @param pathway_names Character vector of pathway names
/// @param gsea_param Weighting parameter (default 1.0)
/// @param n_anchors Number of anchor samples (default: 200 for gamma, 500 for gpd)
/// @param min_size Minimum pathway size (default 1)
/// @param max_size Maximum pathway size (default Inf)
/// @param score_type Score type: "std", "pos", or "neg" (default "std")
/// @param method Method for p-value estimation: "gamma" or "empirical_gpd" (default "gamma")
/// @return Data frame with GSEA results
/// @export
#[extendr]
fn run_gsea_rust_with_method(
    stats: &[f64],
    pathway_indices: List,
    pathway_names: Vec<String>,
    gsea_param: f64,
    n_anchors: i32,
    min_size: i32,
    max_size: i32,
    score_type: &str,
    method: &str,
) -> Robj {
    // Convert pathway indices from R (1-based) to Rust (0-based)
    let pathways: Vec<Vec<usize>> = pathway_indices
        .iter()
        .map(|(_, robj)| {
            if let Some(indices) = robj.as_integer_slice() {
                indices
                    .iter()
                    .map(|&i| (i - 1) as usize)
                    .filter(|&i| i < stats.len())
                    .collect()
            } else {
                Vec::new()
            }
        })
        .collect();

    // Parse method
    let pvalue_method = match method.to_lowercase().as_str() {
        "empirical_gpd" | "gpd" | "empirical" => parallel::PvalueMethod::EmpiricalGpd,
        _ => parallel::PvalueMethod::Gamma, // Default to Gamma
    };

    // Filter pathways by size
    let min_sz = min_size as usize;
    let max_sz = if max_size < 0 { usize::MAX } else { max_size as usize };

    let results = parallel::run_gsea_parallel_with_method(
        stats,
        &pathways,
        &pathway_names,
        gsea_param,
        n_anchors as usize,
        min_sz,
        max_sz,
        score_type,
        pvalue_method,
    );

    // Build output vectors
    let n = results.len();
    let mut out_pathway: Vec<String> = Vec::with_capacity(n);
    let mut out_pval: Vec<f64> = Vec::with_capacity(n);
    let mut out_es: Vec<f64> = Vec::with_capacity(n);
    let mut out_nes: Vec<f64> = Vec::with_capacity(n);
    let mut out_size: Vec<i32> = Vec::with_capacity(n);
    let mut out_leading_edge: Vec<Robj> = Vec::with_capacity(n);

    for result in results {
        out_pathway.push(result.pathway);
        out_pval.push(result.pval);
        out_es.push(result.es);
        out_nes.push(result.nes);
        out_size.push(result.size as i32);
        // Convert leading edge back to 1-based
        let le: Vec<i32> = result.leading_edge.iter().map(|&i| (i + 1) as i32).collect();
        out_leading_edge.push(le.into_robj());
    }

    // Return as a list - R will handle data frame conversion
    list!(
        pathway = out_pathway,
        pval = out_pval,
        ES = out_es,
        NES = out_nes,
        size = out_size,
        leadingEdge = List::from_values(out_leading_edge)
    ).into()
}

/// Run GSEA analysis with adaptive sampling
/// @param stats Named numeric vector of gene statistics
/// @param pathway_indices List of integer vectors, each containing gene indices for a pathway
/// @param pathway_names Character vector of pathway names
/// @param gsea_param Weighting parameter (default 1.0)
/// @param min_size Minimum pathway size (default 1)
/// @param max_size Maximum pathway size (default Inf)
/// @param score_type Score type: "std", "pos", or "neg" (default "std")
/// @return List with GSEA results and total_samples_used
/// @export
#[extendr]
fn run_gsea_rust_adaptive(
    stats: &[f64],
    pathway_indices: List,
    pathway_names: Vec<String>,
    gsea_param: f64,
    min_size: i32,
    max_size: i32,
    score_type: &str,
) -> Robj {
    // Convert pathway indices from R (1-based) to Rust (0-based)
    let pathways: Vec<Vec<usize>> = pathway_indices
        .iter()
        .map(|(_, robj)| {
            if let Some(indices) = robj.as_integer_slice() {
                indices
                    .iter()
                    .map(|&i| (i - 1) as usize)
                    .filter(|&i| i < stats.len())
                    .collect()
            } else {
                Vec::new()
            }
        })
        .collect();

    // Filter pathways by size
    let min_sz = min_size as usize;
    let max_sz = if max_size < 0 { usize::MAX } else { max_size as usize };

    let (results, total_samples) = parallel::run_gsea_parallel_adaptive(
        stats,
        &pathways,
        &pathway_names,
        gsea_param,
        min_sz,
        max_sz,
        score_type,
    );

    // Build output vectors
    let n = results.len();
    let mut out_pathway: Vec<String> = Vec::with_capacity(n);
    let mut out_pval: Vec<f64> = Vec::with_capacity(n);
    let mut out_es: Vec<f64> = Vec::with_capacity(n);
    let mut out_nes: Vec<f64> = Vec::with_capacity(n);
    let mut out_size: Vec<i32> = Vec::with_capacity(n);
    let mut out_leading_edge: Vec<Robj> = Vec::with_capacity(n);

    for result in results {
        out_pathway.push(result.pathway);
        out_pval.push(result.pval);
        out_es.push(result.es);
        out_nes.push(result.nes);
        out_size.push(result.size as i32);
        let le: Vec<i32> = result.leading_edge.iter().map(|&i| (i + 1) as i32).collect();
        out_leading_edge.push(le.into_robj());
    }

    // Return as a list with total_samples
    list!(
        pathway = out_pathway,
        pval = out_pval,
        ES = out_es,
        NES = out_nes,
        size = out_size,
        leadingEdge = List::from_values(out_leading_edge),
        total_samples = total_samples as i32
    ).into()
}

/// Run GSEA with two-pass coarse-to-fine strategy
/// @param stats Named numeric vector of gene statistics
/// @param pathway_indices List of integer vectors
/// @param pathway_names Character vector of pathway names
/// @param gsea_param Weighting parameter (default 1.0)
/// @param min_size Minimum pathway size (default 1)
/// @param max_size Maximum pathway size (default Inf)
/// @param score_type Score type: "std", "pos", or "neg" (default "std")
/// @param coarse_anchors Anchors for coarse pass (default 50)
/// @param fine_anchors Anchors for fine pass (default 200)
/// @param pvalue_threshold P-value threshold for fine pass (default 0.2)
/// @return List with GSEA results and pass statistics
/// @export
#[extendr]
fn run_gsea_rust_two_pass(
    stats: &[f64],
    pathway_indices: List,
    pathway_names: Vec<String>,
    gsea_param: f64,
    min_size: i32,
    max_size: i32,
    score_type: &str,
    coarse_anchors: i32,
    fine_anchors: i32,
    pvalue_threshold: f64,
) -> Robj {
    // Convert pathway indices
    let pathways: Vec<Vec<usize>> = pathway_indices
        .iter()
        .map(|(_, robj)| {
            if let Some(indices) = robj.as_integer_slice() {
                indices
                    .iter()
                    .map(|&i| (i - 1) as usize)
                    .filter(|&i| i < stats.len())
                    .collect()
            } else {
                Vec::new()
            }
        })
        .collect();

    let min_sz = min_size as usize;
    let max_sz = if max_size < 0 { usize::MAX } else { max_size as usize };

    let config = parallel::TwoPassConfig {
        coarse_anchors: coarse_anchors as usize,
        fine_anchors: fine_anchors as usize,
        pvalue_threshold,
    };

    let (results, coarse_count, fine_count) = parallel::run_gsea_parallel_two_pass(
        stats,
        &pathways,
        &pathway_names,
        gsea_param,
        min_sz,
        max_sz,
        score_type,
        &config,
    );

    // Build output
    let n = results.len();
    let mut out_pathway: Vec<String> = Vec::with_capacity(n);
    let mut out_pval: Vec<f64> = Vec::with_capacity(n);
    let mut out_es: Vec<f64> = Vec::with_capacity(n);
    let mut out_nes: Vec<f64> = Vec::with_capacity(n);
    let mut out_size: Vec<i32> = Vec::with_capacity(n);
    let mut out_leading_edge: Vec<Robj> = Vec::with_capacity(n);

    for result in results {
        out_pathway.push(result.pathway);
        out_pval.push(result.pval);
        out_es.push(result.es);
        out_nes.push(result.nes);
        out_size.push(result.size as i32);
        let le: Vec<i32> = result.leading_edge.iter().map(|&i| (i + 1) as i32).collect();
        out_leading_edge.push(le.into_robj());
    }

    list!(
        pathway = out_pathway,
        pval = out_pval,
        ES = out_es,
        NES = out_nes,
        size = out_size,
        leadingEdge = List::from_values(out_leading_edge),
        coarse_filtered = coarse_count as i32,
        fine_calculated = fine_count as i32
    ).into()
}

/// Run GSEA with batch optimization (shared samples + vectorized ES)
/// @param stats Named numeric vector of gene statistics
/// @param pathway_indices List of integer vectors
/// @param pathway_names Character vector of pathway names
/// @param gsea_param Weighting parameter (default 1.0)
/// @param n_anchors Number of anchor samples (default 200)
/// @param min_size Minimum pathway size (default 1)
/// @param max_size Maximum pathway size (default Inf)
/// @param score_type Score type: "std", "pos", or "neg" (default "std")
/// @return Data frame with GSEA results
/// @export
#[extendr]
fn run_gsea_rust_batch(
    stats: &[f64],
    pathway_indices: List,
    pathway_names: Vec<String>,
    gsea_param: f64,
    n_anchors: i32,
    min_size: i32,
    max_size: i32,
    score_type: &str,
) -> Robj {
    // Convert pathway indices from R (1-based) to Rust (0-based)
    let pathways: Vec<Vec<usize>> = pathway_indices
        .iter()
        .map(|(_, robj)| {
            if let Some(indices) = robj.as_integer_slice() {
                indices
                    .iter()
                    .map(|&i| (i - 1) as usize)
                    .filter(|&i| i < stats.len())
                    .collect()
            } else {
                Vec::new()
            }
        })
        .collect();

    let min_sz = min_size as usize;
    let max_sz = if max_size < 0 { usize::MAX } else { max_size as usize };

    let results = parallel::run_gsea_parallel_batch(
        stats,
        &pathways,
        &pathway_names,
        gsea_param,
        n_anchors as usize,
        min_sz,
        max_sz,
        score_type,
    );

    // Build output vectors
    let n = results.len();
    let mut out_pathway: Vec<String> = Vec::with_capacity(n);
    let mut out_pval: Vec<f64> = Vec::with_capacity(n);
    let mut out_es: Vec<f64> = Vec::with_capacity(n);
    let mut out_nes: Vec<f64> = Vec::with_capacity(n);
    let mut out_size: Vec<i32> = Vec::with_capacity(n);
    let mut out_leading_edge: Vec<Robj> = Vec::with_capacity(n);

    for result in results {
        out_pathway.push(result.pathway);
        out_pval.push(result.pval);
        out_es.push(result.es);
        out_nes.push(result.nes);
        out_size.push(result.size as i32);
        let le: Vec<i32> = result.leading_edge.iter().map(|&i| (i + 1) as i32).collect();
        out_leading_edge.push(le.into_robj());
    }

    list!(
        pathway = out_pathway,
        pval = out_pval,
        ES = out_es,
        NES = out_nes,
        size = out_size,
        leadingEdge = List::from_values(out_leading_edge)
    ).into()
}

// Macro to generate exports
extendr_module! {
    mod superfastgsea;
    fn rust_hello;
    fn calc_es_rust;
    fn run_gsea_rust;
    fn run_gsea_rust_with_method;
    fn run_gsea_rust_adaptive;
    fn run_gsea_rust_two_pass;
    fn run_gsea_rust_batch;
}
