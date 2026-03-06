#' Super Fast Gene Set Enrichment Analysis
#'
#' Performs fast GSEA using a Rust backend with distribution approximation.
#' API is compatible with fgsea::fgsea().
#'
#' @param pathways List of gene sets (character vectors of gene names)
#' @param stats Named numeric vector of gene-level statistics (e.g., log2FC or t-statistics)
#' @param minSize Minimum gene set size (default: 1)
#' @param maxSize Maximum gene set size (default: Inf)
#' @param scoreType Score type: "std" (default, two-tailed), "pos" (positive only), "neg" (negative only)
#' @param eps Minimum p-value threshold (default: 1e-10)
#' @param nproc Number of processors to use (0 = auto-detect, default: 0)
#' @param gseaParam GSEA weighting parameter (default: 1)
#' @param nAnchors Number of anchor samples for distribution fitting.
#'   For "gamma" method: default 200, minimum 200.
#'   For "empirical_gpd" method: default 500, minimum 500 (recommended for accuracy).
#' @param method Method for p-value estimation:
#'   \itemize{
#'     \item "gamma": Gamma distribution approximation (fast, default)
#'     \item "gamma_adaptive": Gamma with adaptive sampling (stops when converged)
#'     \item "gamma_two_pass": Two-pass coarse-to-fine (fastest for large datasets)
#'     \item "empirical_gpd": Empirical CDF + GPD for tails (more accurate, uses EVT)
#'     \item "hybrid": Use fgsea for small pathways, Gamma for large (highest accuracy, requires fgsea)
#'   }
#' @param pAdjustMethod Method for multiple testing correction:
#'   \itemize{
#'     \item "BH": Benjamini-Hochberg FDR control (default, standard in bioinformatics)
#'     \item "sidak": Sidak FWER control (more conservative, like blitzGSEA)
#'     \item "bonferroni": Bonferroni FWER control (most conservative)
#'     \item "holm", "hochberg", "hommel", "BY", "none": Other methods from p.adjust()
#'   }
#' @param hybridSize Size threshold for hybrid method. Pathways with size < hybridSize
#'   use fgsea (more accurate for small sets), others use Gamma (faster). Default: 50.
#' @param hybridPvalueThreshold Not used in current implementation (kept for compatibility).
#' @param hybridNperm Number of permutations for fgsea in hybrid method.
#'   Default: 1000. Higher values = more accurate but slower.
#'
#' @return A data.frame with columns:
#'   \item{pathway}{Pathway name}
#'   \item{pval}{Nominal p-value}
#'   \item{padj}{Adjusted p-value (method specified by pAdjustMethod)}
#'   \item{ES}{Enrichment Score}
#'   \item{NES}{Normalized Enrichment Score}
#'   \item{size}{Number of genes in pathway after filtering}
#'   \item{leadingEdge}{List of leading edge gene names}
#'   If pAdjustMethod includes multiple methods, additional columns (padj_sidak, etc.) are added.
#'
#' @details
#' Two methods are available for p-value estimation:
#'
#' \strong{Gamma method (default):}
#' Fast parametric approach that assumes enrichment scores follow a Gamma distribution.
#' Uses 200 anchor samples per size bin. Best for speed when distributional assumption holds.
#'
#' \strong{Empirical+GPD method:}
#' Semi-parametric approach using empirical CDF for the center of the distribution
#' and Generalized Pareto Distribution (GPD) for the tails. Based on Extreme Value Theory.
#' Uses 500+ anchor samples for reliable estimation. More accurate for extreme p-values
#' and when the Gamma assumption may not hold.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' set.seed(42)
#' stats <- setNames(rnorm(1000), paste0("gene", 1:1000))
#' stats <- sort(stats, decreasing = TRUE)
#'
#' pathways <- list(
#'   pathway1 = paste0("gene", 1:50),
#'   pathway2 = paste0("gene", 500:550)
#' )
#'
#' # Fast analysis with Gamma method (default)
#' result_gamma <- superfastGSEA(pathways, stats, method = "gamma")
#'
#' # More accurate analysis with Empirical+GPD method
#' result_gpd <- superfastGSEA(pathways, stats, method = "empirical_gpd")
#' }
#'
#' @export
superfastGSEA <- function(
    pathways,
    stats,
    minSize = 1,
    maxSize = Inf,
    scoreType = c("std", "pos", "neg"),
    eps = 1e-10,
    nproc = 0,
    gseaParam = 1,
    nAnchors = NULL,
    method = c("gamma", "gamma_batch", "gamma_adaptive", "gamma_two_pass", "empirical_gpd", "hybrid"),
    pAdjustMethod = c("BH", "sidak", "bonferroni", "holm", "hochberg", "hommel", "BY", "none"),
    hybridSize = 50,
    hybridPvalueThreshold = 0.1,
    hybridNperm = 1000
) {
    # Validate inputs
    scoreType <- match.arg(scoreType)
    method <- match.arg(method)
    pAdjustMethod <- match.arg(pAdjustMethod)

    # Handle hybrid method separately
    if (method == "hybrid") {
        return(.runHybridGSEA(pathways, stats, minSize, maxSize, scoreType,
                             eps, gseaParam, nAnchors, hybridSize,
                             pvalueThreshold = hybridPvalueThreshold,
                             nPermSimple = hybridNperm,
                             pAdjustMethod = pAdjustMethod))
    }

    # Handle batch optimized method
    if (method == "gamma_batch") {
        return(.runBatchGSEA(pathways, stats, minSize, maxSize, scoreType,
                             eps, gseaParam, nAnchors, pAdjustMethod))
    }

    # Handle adaptive sampling method
    if (method == "gamma_adaptive") {
        return(.runAdaptiveGSEA(pathways, stats, minSize, maxSize, scoreType,
                                eps, gseaParam, pAdjustMethod))
    }

    # Handle two-pass method
    if (method == "gamma_two_pass") {
        return(.runTwoPassGSEA(pathways, stats, minSize, maxSize, scoreType,
                               eps, gseaParam, pAdjustMethod = pAdjustMethod))
    }

    if (!is.list(pathways)) {
        stop("'pathways' must be a list of character vectors")
    }

    if (!is.numeric(stats) || is.null(names(stats))) {
        stop("'stats' must be a named numeric vector")
    }

    # Set default nAnchors based on method
    if (is.null(nAnchors)) {
        nAnchors <- switch(method,
            "gamma" = 200L,
            "empirical_gpd" = 500L,
            200L  # default
        )
    }

    # Prepare shared inputs
    inputs <- .prepareGSEAInputs(pathways, stats)
    maxSize_int <- if (is.infinite(maxSize)) -1L else as.integer(maxSize)

    # Call Rust backend with method selection
    result <- run_gsea_rust_with_method(
        stats = inputs$stats_values,
        pathway_indices = inputs$pathway_indices,
        pathway_names = inputs$pathway_names,
        gsea_param = as.double(gseaParam),
        n_anchors = as.integer(nAnchors),
        min_size = as.integer(minSize),
        max_size = maxSize_int,
        score_type = scoreType,
        method = method
    )

    .formatGSEAResult(result, inputs$gene_names, eps, pAdjustMethod)
}

#' Test Rust Bridge
#'
#' Simple test function to verify the Rust backend is working.
#'
#' @return A greeting message from Rust
#' @export
test_rust_bridge <- function() {
    rust_hello()
}

#' Calculate Enrichment Score
#'
#' Calculate ES for a single gene set (low-level function).
#'
#' @param stats Named numeric vector of gene statistics (sorted)
#' @param geneSet Character vector of gene names in the set
#' @param gseaParam GSEA weighting parameter (default: 1)
#'
#' @return List with ES and leading edge genes
#' @export
calcES <- function(stats, geneSet, gseaParam = 1) {
    if (!is.numeric(stats) || is.null(names(stats))) {
        stop("'stats' must be a named numeric vector")
    }

    # Sort stats
    stats <- sort(stats, decreasing = TRUE)
    gene_names <- names(stats)
    stats_values <- unname(stats)

    # Get indices of genes in the set
    gene_to_idx <- setNames(seq_along(gene_names), gene_names)
    indices <- as.integer(gene_to_idx[geneSet])
    indices <- indices[!is.na(indices)]

    if (length(indices) == 0) {
        return(list(ES = 0, leadingEdge = character(0)))
    }

    # Call Rust function
    result <- calc_es_rust(stats_values, indices, as.double(gseaParam))

    # Convert leading edge indices to gene names
    result$leadingEdge <- gene_names[result$leadingEdge]

    return(result)
}

#' Internal: Run GSEA with two-pass coarse-to-fine strategy
#' @keywords internal
.runTwoPassGSEA <- function(pathways, stats, minSize, maxSize, scoreType,
                             eps, gseaParam,
                             coarseAnchors = 50L,
                             fineAnchors = 200L,
                             pvalueThreshold = 0.2,
                             pAdjustMethod = "BH") {
    inputs <- .prepareGSEAInputs(pathways, stats)
    maxSize_int <- if (is.infinite(maxSize)) -1L else as.integer(maxSize)

    # Call Rust backend with two-pass
    result <- run_gsea_rust_two_pass(
        stats = inputs$stats_values,
        pathway_indices = inputs$pathway_indices,
        pathway_names = inputs$pathway_names,
        gsea_param = as.double(gseaParam),
        min_size = as.integer(minSize),
        max_size = maxSize_int,
        score_type = scoreType,
        coarse_anchors = as.integer(coarseAnchors),
        fine_anchors = as.integer(fineAnchors),
        pvalue_threshold = as.double(pvalueThreshold)
    )

    .formatGSEAResult(result, inputs$gene_names, eps, pAdjustMethod,
                      extra_fields = c("coarse_filtered", "fine_calculated"))
}

#' Internal: Run GSEA with adaptive sampling
#' @keywords internal
.runAdaptiveGSEA <- function(pathways, stats, minSize, maxSize, scoreType,
                              eps, gseaParam, pAdjustMethod = "BH") {
    inputs <- .prepareGSEAInputs(pathways, stats)
    maxSize_int <- if (is.infinite(maxSize)) -1L else as.integer(maxSize)

    # Call Rust backend with adaptive sampling
    result <- run_gsea_rust_adaptive(
        stats = inputs$stats_values,
        pathway_indices = inputs$pathway_indices,
        pathway_names = inputs$pathway_names,
        gsea_param = as.double(gseaParam),
        min_size = as.integer(minSize),
        max_size = maxSize_int,
        score_type = scoreType
    )

    .formatGSEAResult(result, inputs$gene_names, eps, pAdjustMethod,
                      extra_fields = c("total_samples"))
}

#' Internal: Run hybrid GSEA (fgsea for small pathways, Gamma for large pathways)
#' @keywords internal
#' @param hybridSize Size threshold: pathways < hybridSize use fgsea, >= hybridSize use Gamma
#' @param nPermSimple Number of permutations for fgsea (default 1000)
.runHybridGSEA <- function(pathways, stats, minSize, maxSize, scoreType,
                           eps, gseaParam, nAnchors, hybridSize,
                           pvalueThreshold = 0.1, nPermSimple = 1000,
                           pAdjustMethod = "BH") {
    # Check if fgsea is available
    if (!requireNamespace("fgsea", quietly = TRUE)) {
        stop("Package 'fgsea' is required for hybrid method. Install with: BiocManager::install('fgsea')")
    }

    inputs <- .prepareGSEAInputs(pathways, stats)
    gene_names <- inputs$gene_names
    # Reconstruct named stats for fgsea (requires named vector)
    stats <- setNames(inputs$stats_values, gene_names)

    # Calculate pathway sizes (after filtering for genes in stats)
    pathway_sizes <- sapply(pathways, function(genes) {
        sum(genes %in% gene_names)
    })

    # Split pathways by size
    small_idx <- which(pathway_sizes >= minSize & pathway_sizes < hybridSize &
                       (is.infinite(maxSize) | pathway_sizes <= maxSize))
    large_idx <- which(pathway_sizes >= hybridSize &
                       (is.infinite(maxSize) | pathway_sizes <= maxSize))

    result_list <- list()

    # Process small pathways with fgsea (more accurate for small sets)
    if (length(small_idx) > 0) {
        small_pathways <- pathways[small_idx]
        result_fgsea <- fgsea::fgsea(
            pathways = small_pathways,
            stats = stats,
            minSize = minSize,
            maxSize = hybridSize - 1,
            scoreType = scoreType,
            eps = eps,
            gseaParam = gseaParam,
            nPermSimple = nPermSimple
        )
        result_fgsea <- as.data.frame(result_fgsea)
        result_list$small <- result_fgsea
    }

    # Process large pathways with Gamma (fast and accurate for large sets)
    if (length(large_idx) > 0) {
        large_pathways <- pathways[large_idx]
        result_gamma <- superfastGSEA(
            pathways = large_pathways,
            stats = stats,
            minSize = hybridSize,
            maxSize = maxSize,
            scoreType = scoreType,
            eps = eps,
            gseaParam = gseaParam,
            nAnchors = if (is.null(nAnchors)) 200L else nAnchors,
            method = "gamma"
        )
        result_list$large <- result_gamma
    }

    # Combine results
    if (length(result_list) == 0) {
        return(data.frame(
            pathway = character(0),
            pval = numeric(0),
            padj = numeric(0),
            ES = numeric(0),
            NES = numeric(0),
            size = integer(0),
            leadingEdge = I(list())
        ))
    }

    # Standardize columns and combine
    standardize_result <- function(df) {
        df[, c("pathway", "pval", "ES", "NES", "size", "leadingEdge"), drop = FALSE]
    }

    result_list <- lapply(result_list, standardize_result)
    result <- do.call(rbind, result_list)

    # Recalculate adjusted p-values
    result$padj <- .adjustPvalues(result$pval, pAdjustMethod)
    result$padj <- pmax(result$padj, eps)

    # Sort by p-value
    result <- result[order(result$pval), ]
    rownames(result) <- NULL

    # Reorder columns
    result <- result[, c("pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge")]

    # Add attribute
    attr(result, "pAdjustMethod") <- pAdjustMethod

    return(result)
}

#' Internal: Run GSEA with batch optimization
#' @keywords internal
.runBatchGSEA <- function(pathways, stats, minSize, maxSize, scoreType,
                          eps, gseaParam, nAnchors, pAdjustMethod = "BH") {
    inputs <- .prepareGSEAInputs(pathways, stats)
    maxSize_int <- if (is.infinite(maxSize)) -1L else as.integer(maxSize)

    # Set default nAnchors
    if (is.null(nAnchors)) {
        nAnchors <- 200L
    }

    # Call Rust backend with batch optimization
    result <- run_gsea_rust_batch(
        stats = inputs$stats_values,
        pathway_indices = inputs$pathway_indices,
        pathway_names = inputs$pathway_names,
        gsea_param = as.double(gseaParam),
        n_anchors = as.integer(nAnchors),
        min_size = as.integer(minSize),
        max_size = maxSize_int,
        score_type = scoreType
    )

    .formatGSEAResult(result, inputs$gene_names, eps, pAdjustMethod)
}

#' Internal: Prepare GSEA inputs (shared by all method runners)
#' @keywords internal
.prepareGSEAInputs <- function(pathways, stats) {
    # Remove NA values with warning
    na_count <- sum(is.na(stats))
    if (na_count > 0) {
        warning(sprintf("Removed %d genes with NA statistics", na_count))
        stats <- stats[!is.na(stats)]
    }

    # Sort stats in decreasing order
    stats <- sort(stats, decreasing = TRUE)
    gene_names <- names(stats)
    stats_values <- unname(stats)

    # Convert pathways to indices
    pathway_names <- names(pathways)
    if (is.null(pathway_names)) {
        pathway_names <- paste0("pathway_", seq_along(pathways))
    }

    # Map gene names to indices (1-based for R)
    gene_to_idx <- setNames(seq_along(gene_names), gene_names)
    pathway_indices <- lapply(pathways, function(genes) {
        idx <- gene_to_idx[genes]
        as.integer(idx[!is.na(idx)])
    })

    list(
        gene_names = gene_names,
        stats_values = stats_values,
        pathway_names = pathway_names,
        pathway_indices = pathway_indices
    )
}

#' Internal: Format Rust GSEA result into standard data.frame
#' @keywords internal
.formatGSEAResult <- function(result, gene_names, eps, pAdjustMethod,
                               extra_fields = NULL) {
    # Extract leadingEdge separately as it's a list column
    le <- result$leadingEdge
    result$leadingEdge <- NULL

    # Remove any extra fields from the list before making data.frame
    extra_values <- list()
    if (!is.null(extra_fields)) {
        for (f in extra_fields) {
            extra_values[[f]] <- result[[f]]
            result[[f]] <- NULL
        }
    }

    # Create data frame from other columns
    result <- data.frame(
        pathway = result$pathway,
        pval = result$pval,
        ES = result$ES,
        NES = result$NES,
        size = result$size,
        stringsAsFactors = FALSE
    )
    result$leadingEdge <- le

    # Calculate adjusted p-values
    result$padj <- .adjustPvalues(result$pval, pAdjustMethod)

    # Apply eps threshold
    result$pval <- pmax(result$pval, eps)
    result$padj <- pmax(result$padj, eps)

    # Convert leading edge indices to gene names
    result$leadingEdge <- lapply(result$leadingEdge, function(idx) {
        gene_names[idx]
    })

    # Reorder columns to match fgsea output
    result <- result[, c("pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge")]

    # Sort by p-value
    result <- result[order(result$pval), ]
    rownames(result) <- NULL

    # Attach extra fields as attributes
    for (f in names(extra_values)) {
        attr(result, f) <- extra_values[[f]]
    }
    attr(result, "pAdjustMethod") <- pAdjustMethod

    result
}

#' Internal: Calculate adjusted p-values with multiple methods
#' @param pval Numeric vector of raw p-values
#' @param method Adjustment method: "BH", "sidak", "bonferroni", etc.
#' @return Numeric vector of adjusted p-values
#' @keywords internal
.adjustPvalues <- function(pval, method = "BH") {
    n <- length(pval)

    if (method == "sidak") {
        # Sidak correction: p_adj = 1 - (1 - p)^n
        # More powerful than Bonferroni but still controls FWER
        padj <- 1 - (1 - pval)^n
        # Ensure values are in [0, 1]
        padj <- pmin(pmax(padj, 0), 1)
    } else if (method %in% c("BH", "bonferroni", "holm", "hochberg", "hommel", "BY", "none")) {
        # Use stats::p.adjust for standard methods
        padj <- stats::p.adjust(pval, method = method)
    } else {
        warning(sprintf("Unknown adjustment method '%s', using 'BH'", method))
        padj <- stats::p.adjust(pval, method = "BH")
    }

    return(padj)
}
