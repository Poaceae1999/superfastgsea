# Statistical Validation: superfastGSEA vs fgsea
#
# This script validates that superfastGSEA produces statistically
# equivalent results to fgsea.

library(superfastgsea)

# Check if fgsea is available
if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea package is required for validation. Install with: BiocManager::install('fgsea')")
}
library(fgsea)

#' Calculate concordance metrics between two GSEA results
#'
#' @param result1 First GSEA result data frame
#' @param result2 Second GSEA result data frame
#' @return List of concordance metrics
calculate_concordance <- function(result1, result2) {
    # Merge results by pathway name
    merged <- merge(result1, result2, by = "pathway", suffixes = c(".super", ".fgsea"))

    if (nrow(merged) == 0) {
        return(list(
            nes_cor = NA,
            pval_cor = NA,
            es_cor = NA,
            top50_overlap = NA,
            n_common = 0
        ))
    }

    # NES correlation
    nes_cor <- cor(merged$NES.super, merged$NES.fgsea, method = "spearman", use = "complete.obs")

    # P-value correlation (on -log10 scale)
    pval_cor <- cor(
        -log10(pmax(merged$pval.super, 1e-100)),
        -log10(pmax(merged$pval.fgsea, 1e-100)),
        method = "spearman",
        use = "complete.obs"
    )

    # ES correlation
    es_cor <- cor(merged$ES.super, merged$ES.fgsea, method = "spearman", use = "complete.obs")

    # Top 50 pathway overlap
    top50_super <- head(result1[order(result1$pval), "pathway"], 50)
    top50_fgsea <- head(result2[order(result2$pval), "pathway"], 50)
    top50_overlap <- length(intersect(top50_super, top50_fgsea)) / 50

    # NES difference statistics
    nes_diff <- merged$NES.super - merged$NES.fgsea
    nes_rmsd <- sqrt(mean(nes_diff^2, na.rm = TRUE))
    nes_mad <- median(abs(nes_diff), na.rm = TRUE)

    return(list(
        nes_cor = nes_cor,
        pval_cor = pval_cor,
        es_cor = es_cor,
        top50_overlap = top50_overlap,
        nes_rmsd = nes_rmsd,
        nes_mad = nes_mad,
        n_common = nrow(merged)
    ))
}

#' Run validation against fgsea
#'
#' @param stats Named numeric vector of gene statistics
#' @param pathways List of gene sets
#' @param verbose Print progress
#' @return Concordance metrics
validate_against_fgsea <- function(stats, pathways, verbose = TRUE) {
    if (verbose) cat("Running superfastGSEA...\n")
    result_super <- superfastGSEA(pathways, stats, minSize = 15, maxSize = 500)

    if (verbose) cat("Running fgsea...\n")
    result_fgsea <- fgsea(pathways, stats, minSize = 15, maxSize = 500, nPermSimple = 10000)
    result_fgsea <- as.data.frame(result_fgsea)

    if (verbose) cat("Calculating concordance...\n")
    concordance <- calculate_concordance(result_super, result_fgsea)

    return(concordance)
}

#' Generate synthetic data with known enrichment
#'
#' @param n_genes Total number of genes
#' @param n_pathways Number of pathways
#' @param n_enriched Number of truly enriched pathways
#' @param effect_size Effect size for enriched pathways
#' @return List with stats and pathways
generate_synthetic_data <- function(n_genes = 10000, n_pathways = 100,
                                     n_enriched = 10, effect_size = 2) {
    gene_names <- paste0("gene", 1:n_genes)

    # Generate pathways
    pathways <- lapply(1:n_pathways, function(i) {
        size <- sample(20:200, 1)
        sample(gene_names, size)
    })
    names(pathways) <- paste0("pathway", 1:n_pathways)

    # Generate statistics (null distribution)
    stats <- setNames(rnorm(n_genes), gene_names)

    # Inject signal into first n_enriched pathways
    for (i in 1:n_enriched) {
        pathway_genes <- pathways[[i]]
        idx <- match(pathway_genes, gene_names)
        stats[idx] <- stats[idx] + effect_size
    }

    return(list(
        stats = sort(stats, decreasing = TRUE),
        pathways = pathways,
        enriched = paste0("pathway", 1:n_enriched)
    ))
}

#' Validate power and FDR control
#'
#' @param n_iterations Number of simulation iterations
#' @return Summary statistics
validate_power_fdr <- function(n_iterations = 10) {
    results <- lapply(1:n_iterations, function(i) {
        cat(sprintf("Iteration %d/%d\n", i, n_iterations))

        # Generate data
        data <- generate_synthetic_data(n_enriched = 10, effect_size = 2.5)

        # Run superfastGSEA
        result <- superfastGSEA(data$pathways, data$stats, minSize = 15)

        # Calculate metrics
        significant <- result$pathway[result$padj < 0.25]
        true_positives <- length(intersect(significant, data$enriched))
        false_positives <- length(setdiff(significant, data$enriched))

        # Power: proportion of enriched pathways detected
        power <- true_positives / length(data$enriched)

        # FDR: proportion of false positives among significant
        actual_fdr <- if (length(significant) > 0) {
            false_positives / length(significant)
        } else {
            0
        }

        return(list(power = power, actual_fdr = actual_fdr))
    })

    # Summarize
    powers <- sapply(results, function(r) r$power)
    fdrs <- sapply(results, function(r) r$actual_fdr)

    cat("\n=== Power and FDR Summary ===\n")
    cat(sprintf("Mean Power: %.2f (SD: %.2f)\n", mean(powers), sd(powers)))
    cat(sprintf("Mean Actual FDR: %.2f (SD: %.2f)\n", mean(fdrs), sd(fdrs)))
    cat(sprintf("FDR Control: %s (target: 0.25)\n",
                if (mean(fdrs) <= 0.25) "PASS" else "FAIL"))

    return(list(powers = powers, fdrs = fdrs))
}

# Main execution
if (interactive()) {
    cat("=== superfastGSEA Statistical Validation ===\n\n")

    # Test 1: Concordance with fgsea on real-like data
    cat("--- Test 1: Concordance with fgsea ---\n")
    set.seed(42)
    data <- generate_synthetic_data(n_genes = 10000, n_pathways = 200)
    concordance <- validate_against_fgsea(data$stats, data$pathways)

    cat("\nConcordance Metrics:\n")
    cat(sprintf("  NES Spearman correlation: %.4f (target: > 0.99)\n", concordance$nes_cor))
    cat(sprintf("  P-value Spearman correlation: %.4f (target: > 0.98)\n", concordance$pval_cor))
    cat(sprintf("  ES Spearman correlation: %.4f\n", concordance$es_cor))
    cat(sprintf("  Top 50 pathway overlap: %.2f (target: > 0.95)\n", concordance$top50_overlap))
    cat(sprintf("  NES median absolute difference: %.4f (target: < 0.05)\n", concordance$nes_mad))

    # Validation status
    pass_nes <- concordance$nes_cor > 0.99
    pass_pval <- concordance$pval_cor > 0.98
    pass_top50 <- concordance$top50_overlap > 0.95
    pass_nes_mad <- concordance$nes_mad < 0.05

    cat(sprintf("\nOverall: %s\n",
                if (all(c(pass_nes, pass_pval, pass_top50))) "PASS" else "NEEDS REVIEW"))

    # Test 2: Power and FDR (uncomment to run full test)
    # cat("\n--- Test 2: Power and FDR Control ---\n")
    # power_results <- validate_power_fdr(n_iterations = 10)
}
