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

// ---------------------------------------------------------------------------
// Shared FFI helpers
// ---------------------------------------------------------------------------

/// Convert R 1-based pathway index list to Rust 0-based Vec<Vec<usize>>.
fn convert_pathway_indices(pathway_indices: &List, stats_len: usize) -> Vec<Vec<usize>> {
    pathway_indices
        .iter()
        .map(|(_, robj)| {
            if let Some(indices) = robj.as_integer_slice() {
                indices
                    .iter()
                    .map(|&i| (i - 1) as usize)
                    .filter(|&i| i < stats_len)
                    .collect()
            } else {
                Vec::new()
            }
        })
        .collect()
}

/// Parse R min/max size into Rust usizes (max_size < 0 → usize::MAX).
fn parse_size_bounds(min_size: i32, max_size: i32) -> (usize, usize) {
    let min_sz = min_size as usize;
    let max_sz = if max_size < 0 { usize::MAX } else { max_size as usize };
    (min_sz, max_sz)
}

/// Pre-built output vectors ready for composing into an R list.
struct GseaOutputVecs {
    pathway: Vec<String>,
    pval: Vec<f64>,
    es: Vec<f64>,
    nes: Vec<f64>,
    size: Vec<i32>,
    leading_edge: Vec<Robj>,
}

/// Collect results into output vectors, converting leading-edge to 1-based.
fn collect_gsea_results(results: Vec<gsea_core::GseaResult>) -> GseaOutputVecs {
    let n = results.len();
    let mut out = GseaOutputVecs {
        pathway: Vec::with_capacity(n),
        pval: Vec::with_capacity(n),
        es: Vec::with_capacity(n),
        nes: Vec::with_capacity(n),
        size: Vec::with_capacity(n),
        leading_edge: Vec::with_capacity(n),
    };

    for r in results {
        out.pathway.push(r.pathway);
        out.pval.push(r.pval);
        out.es.push(r.es);
        out.nes.push(r.nes);
        out.size.push(r.size as i32);
        let le: Vec<i32> = r.leading_edge.iter().map(|&i| (i + 1) as i32).collect();
        out.leading_edge.push(le.into_robj());
    }

    out
}

/// Build the standard R list (pathway, pval, ES, NES, size, leadingEdge).
fn build_gsea_list(out: &GseaOutputVecs) -> List {
    list!(
        pathway = &out.pathway,
        pval = &out.pval,
        ES = &out.es,
        NES = &out.nes,
        size = &out.size,
        leadingEdge = List::from_values(&out.leading_edge)
    )
}

// ---------------------------------------------------------------------------
// Exported FFI functions
// ---------------------------------------------------------------------------

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
    let indices: Vec<usize> = gene_set_indices
        .iter()
        .map(|&i| (i - 1) as usize)
        .filter(|&i| i < stats.len())
        .collect();

    let (es, leading_edge) = gsea_core::calc_enrichment_score_sparse(stats, &indices, gsea_param);
    let leading_edge_r: Vec<i32> = leading_edge.iter().map(|&i| (i + 1) as i32).collect();

    list!(
        ES = es,
        leadingEdge = leading_edge_r
    )
}

/// Run GSEA analysis on multiple pathways in parallel
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
    let pathways = convert_pathway_indices(&pathway_indices, stats.len());
    let (min_sz, max_sz) = parse_size_bounds(min_size, max_size);

    let results = parallel::run_gsea_parallel(
        stats, &pathways, &pathway_names, gsea_param,
        n_anchors as usize, min_sz, max_sz, score_type,
    );

    let out = collect_gsea_results(results);
    build_gsea_list(&out).into()
}

/// Run GSEA analysis with method selection
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
    let pathways = convert_pathway_indices(&pathway_indices, stats.len());
    let (min_sz, max_sz) = parse_size_bounds(min_size, max_size);

    let pvalue_method = match method.to_lowercase().as_str() {
        "empirical_gpd" | "gpd" | "empirical" => parallel::PvalueMethod::EmpiricalGpd,
        _ => parallel::PvalueMethod::Gamma,
    };

    let results = parallel::run_gsea_parallel_with_method(
        stats, &pathways, &pathway_names, gsea_param,
        n_anchors as usize, min_sz, max_sz, score_type, pvalue_method,
    );

    let out = collect_gsea_results(results);
    build_gsea_list(&out).into()
}

/// Run GSEA analysis with adaptive sampling
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
    let pathways = convert_pathway_indices(&pathway_indices, stats.len());
    let (min_sz, max_sz) = parse_size_bounds(min_size, max_size);

    let (results, total_samples) = parallel::run_gsea_parallel_adaptive(
        stats, &pathways, &pathway_names, gsea_param,
        min_sz, max_sz, score_type,
    );

    let out = collect_gsea_results(results);
    list!(
        pathway = &out.pathway,
        pval = &out.pval,
        ES = &out.es,
        NES = &out.nes,
        size = &out.size,
        leadingEdge = List::from_values(&out.leading_edge),
        total_samples = total_samples as i32
    ).into()
}

/// Run GSEA with two-pass coarse-to-fine strategy
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
    let pathways = convert_pathway_indices(&pathway_indices, stats.len());
    let (min_sz, max_sz) = parse_size_bounds(min_size, max_size);

    let config = parallel::TwoPassConfig {
        coarse_anchors: coarse_anchors as usize,
        fine_anchors: fine_anchors as usize,
        pvalue_threshold,
    };

    let (results, coarse_count, fine_count) = parallel::run_gsea_parallel_two_pass(
        stats, &pathways, &pathway_names, gsea_param,
        min_sz, max_sz, score_type, &config,
    );

    let out = collect_gsea_results(results);
    list!(
        pathway = &out.pathway,
        pval = &out.pval,
        ES = &out.es,
        NES = &out.nes,
        size = &out.size,
        leadingEdge = List::from_values(&out.leading_edge),
        coarse_filtered = coarse_count as i32,
        fine_calculated = fine_count as i32
    ).into()
}

/// Run GSEA with batch optimization (shared samples + vectorized ES)
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
    let pathways = convert_pathway_indices(&pathway_indices, stats.len());
    let (min_sz, max_sz) = parse_size_bounds(min_size, max_size);

    let results = parallel::run_gsea_parallel_batch(
        stats, &pathways, &pathway_names, gsea_param,
        n_anchors as usize, min_sz, max_sz, score_type,
    );

    let out = collect_gsea_results(results);
    build_gsea_list(&out).into()
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
