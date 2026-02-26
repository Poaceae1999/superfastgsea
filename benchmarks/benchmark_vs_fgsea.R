# Benchmark: superfastGSEA vs fgsea
#
# This script compares the performance of superfastGSEA against fgsea
# across various data sizes and configurations.

library(superfastgsea)
library(microbenchmark)

# Check if fgsea is available
if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea package is required for benchmarking. Install with: BiocManager::install('fgsea')")
}
library(fgsea)

#' Generate random pathways
#'
#' @param n_pathways Number of pathways to generate
#' @param n_genes Total number of genes in the universe
#' @param size_range Vector of two integers: min and max pathway size
#' @return List of pathways (character vectors)
generate_random_pathways <- function(n_pathways, n_genes, size_range = c(15, 500)) {
    gene_names <- paste0("gene", 1:n_genes)
    pathways <- lapply(1:n_pathways, function(i) {
        size <- sample(size_range[1]:size_range[2], 1)
        sample(gene_names, size)
    })
    names(pathways) <- paste0("pathway", 1:n_pathways)
    return(pathways)
}

#' Run a single benchmark comparison
#'
#' @param n_genes Number of genes
#' @param n_pathways Number of pathways
#' @param size_range Pathway size range
#' @param times Number of benchmark repetitions
#' @return Benchmark results
run_benchmark <- function(n_genes, n_pathways, size_range = c(15, 200), times = 5) {
    cat(sprintf("\nBenchmarking: %d genes, %d pathways, size %d-%d\n",
                n_genes, n_pathways, size_range[1], size_range[2]))

    # Generate test data
    set.seed(42)
    stats <- setNames(rnorm(n_genes), paste0("gene", 1:n_genes))
    pathways <- generate_random_pathways(n_pathways, n_genes, size_range)

    # Run benchmark
    results <- microbenchmark(
        superfastGSEA = {
            superfastGSEA(pathways, stats, minSize = size_range[1])
        },
        fgsea = {
            fgsea(pathways, stats, minSize = size_range[1], nPermSimple = 1000)
        },
        times = times
    )

    # Calculate speedup
    summary_df <- summary(results)
    fgsea_time <- summary_df$median[summary_df$expr == "fgsea"]
    super_time <- summary_df$median[summary_df$expr == "superfastGSEA"]
    speedup <- fgsea_time / super_time

    cat(sprintf("  fgsea median: %.2f ms\n", fgsea_time))
    cat(sprintf("  superfastGSEA median: %.2f ms\n", super_time))
    cat(sprintf("  Speedup: %.1fx\n", speedup))

    return(list(
        results = results,
        speedup = speedup,
        config = list(n_genes = n_genes, n_pathways = n_pathways, size_range = size_range)
    ))
}

#' Run full benchmark matrix
#'
#' @return Data frame with all benchmark results
run_benchmark_matrix <- function() {
    configs <- list(
        # Small scale
        list(n_genes = 5000, n_pathways = 50, size_range = c(15, 100)),
        # Medium scale
        list(n_genes = 10000, n_pathways = 200, size_range = c(15, 200)),
        # Large scale
        list(n_genes = 20000, n_pathways = 500, size_range = c(15, 300)),
        # Very large (stress test)
        list(n_genes = 20000, n_pathways = 1000, size_range = c(15, 500))
    )

    results <- lapply(configs, function(cfg) {
        tryCatch({
            run_benchmark(cfg$n_genes, cfg$n_pathways, cfg$size_range)
        }, error = function(e) {
            cat(sprintf("Error in config: %s\n", e$message))
            NULL
        })
    })

    # Compile results
    summary_df <- data.frame(
        n_genes = sapply(results, function(r) if (!is.null(r)) r$config$n_genes else NA),
        n_pathways = sapply(results, function(r) if (!is.null(r)) r$config$n_pathways else NA),
        speedup = sapply(results, function(r) if (!is.null(r)) r$speedup else NA)
    )

    return(summary_df)
}

# Main execution
if (interactive()) {
    cat("=== superfastGSEA Benchmark Suite ===\n")

    # Quick test
    cat("\n--- Quick Test (50 pathways) ---\n")
    quick_result <- run_benchmark(5000, 50, c(15, 100), times = 3)

    # Full matrix (uncomment to run)
    # cat("\n--- Full Benchmark Matrix ---\n")
    # full_results <- run_benchmark_matrix()
    # print(full_results)
}
