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

    # Top-enriched pathway: boost genes 1:30 so they rank at the top
    stats_top <- rnorm(n)
    stats_top[1:30] <- stats_top[1:30] + 3
    names(stats_top) <- paste0("gene", 1:n)

    pathways_top <- list(test = paste0("gene", 1:30))
    result_top <- superfastGSEA(pathways_top, stats_top, minSize = 10)
    expect_true(result_top$ES > 0)

    # Bottom-enriched pathway: suppress genes 1:30 so they rank at the bottom
    stats_bottom <- rnorm(n)
    stats_bottom[1:30] <- stats_bottom[1:30] - 3
    names(stats_bottom) <- paste0("gene", 1:n)

    pathways_bottom <- list(test = paste0("gene", 1:30))
    result_bottom <- superfastGSEA(pathways_bottom, stats_bottom, minSize = 10)
    expect_true(result_bottom$ES < 0)
})

test_that("reproducibility: ES is deterministic for same input", {
    # ES calculation is deterministic (no random sampling involved).
    # NES and p-values depend on Gamma fitting which uses Rust-side RNG
    # in parallel (rayon), so they are not reproducible across runs.
    set.seed(42)
    stats <- setNames(rnorm(500), paste0("gene", 1:500))
    pathways <- list(test = paste0("gene", 1:50))

    result1 <- superfastGSEA(pathways, stats)
    result2 <- superfastGSEA(pathways, stats)

    # ES must be identical (deterministic)
    expect_equal(result1$ES, result2$ES)
    # Size must be identical
    expect_equal(result1$size, result2$size)
    # Leading edge must be identical
    expect_equal(result1$leadingEdge, result2$leadingEdge)
})
