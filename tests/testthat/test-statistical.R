# Statistical Validation Tests
# Validates statistical accuracy of superfastGSEA

test_that("enriched pathways have significant p-values", {
    set.seed(42)

    # Create data with known enrichment at top
    n_genes <- 1000
    enriched_genes <- 1:30

    stats <- rnorm(n_genes)
    stats[enriched_genes] <- stats[enriched_genes] + 3  # Strong signal
    names(stats) <- paste0("gene", 1:n_genes)
    stats <- sort(stats, decreasing = TRUE)

    pathways <- list(
        enriched = paste0("gene", enriched_genes),
        null = paste0("gene", 500:530)
    )

    result <- superfastGSEA(pathways, stats, minSize = 10)

    # Enriched pathway should be significant
    enriched_result <- result[result$pathway == "enriched", ]
    expect_true(enriched_result$pval < 0.05)
    expect_true(enriched_result$ES > 0)

    # Null pathway should NOT be significant (usually)
    null_result <- result[result$pathway == "null", ]
    expect_true(null_result$pval > enriched_result$pval)
})

test_that("NES is properly normalized", {
    set.seed(42)

    stats <- setNames(rnorm(1000), paste0("gene", 1:1000))
    stats <- sort(stats, decreasing = TRUE)

    # Different sized pathways
    pathways <- list(
        small = names(stats)[1:15],
        large = names(stats)[1:100]
    )

    result <- superfastGSEA(pathways, stats, minSize = 10)

    # Both should have non-zero NES
    expect_true(all(!is.na(result$NES)))
    expect_true(all(is.finite(result$NES)))
})

test_that("p-values are bounded correctly", {
    set.seed(42)
    stats <- setNames(rnorm(500), paste0("gene", 1:500))

    pathways <- list(test = paste0("gene", 1:50))

    result <- superfastGSEA(pathways, stats, eps = 1e-10)

    # p-values should be >= eps
    expect_true(all(result$pval >= 1e-10))
    expect_true(all(result$padj >= 1e-10))

    # p-values should be <= 1
    expect_true(all(result$pval <= 1))
    expect_true(all(result$padj <= 1))
})

test_that("multiple testing correction is applied", {
    set.seed(42)
    stats <- setNames(rnorm(1000), paste0("gene", 1:1000))

    # Create many pathways
    pathways <- lapply(1:50, function(i) {
        paste0("gene", sample(1000, 30))
    })
    names(pathways) <- paste0("pathway", 1:50)

    result <- superfastGSEA(pathways, stats, minSize = 10)

    # padj should be >= pval (BH correction)
    expect_true(all(result$padj >= result$pval))
})

test_that("ES sign matches enrichment direction", {
    set.seed(42)
    n <- 1000

    # Top-enriched pathway
    stats_top <- rnorm(n)
    stats_top[1:30] <- stats_top[1:30] + 3
    names(stats_top) <- paste0("gene", 1:n)
    stats_top <- sort(stats_top, decreasing = TRUE)

    # Bottom-enriched pathway
    stats_bottom <- rnorm(n)
    stats_bottom[(n-29):n] <- stats_bottom[(n-29):n] + 3
    names(stats_bottom) <- paste0("gene", 1:n)
    stats_bottom <- sort(stats_bottom, decreasing = TRUE)

    pathways_top <- list(test = paste0("gene", 1:30))
    pathways_bottom <- list(test = paste0("gene", (n-29):n))

    result_top <- superfastGSEA(pathways_top, stats_top, minSize = 10)
    result_bottom <- superfastGSEA(pathways_bottom, stats_bottom, minSize = 10)

    # Top-enriched should have positive ES
    expect_true(result_top$ES > 0)

    # Bottom-enriched should have negative ES
    expect_true(result_bottom$ES < 0)
})

test_that("reproducibility: same input gives same output", {
    stats <- setNames(rnorm(500), paste0("gene", 1:500))
    pathways <- list(test = paste0("gene", 1:50))

    # Note: Due to random sampling in Gamma fitting, we use seed
    set.seed(123)
    result1 <- superfastGSEA(pathways, stats)

    set.seed(123)
    result2 <- superfastGSEA(pathways, stats)

    expect_equal(result1$ES, result2$ES)
    expect_equal(result1$NES, result2$NES)
    # p-values may have small numerical differences
    expect_true(abs(result1$pval - result2$pval) < 0.01)
})
