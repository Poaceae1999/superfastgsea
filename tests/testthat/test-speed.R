# Speed Benchmark Tests
# Basic speed tests (detailed benchmarks in benchmarks/ directory)

test_that("handles large number of pathways efficiently", {
    skip_on_cran()

    set.seed(42)
    n_genes <- 5000
    n_pathways <- 100

    stats <- setNames(rnorm(n_genes), paste0("gene", 1:n_genes))

    # Generate random pathways
    pathways <- lapply(1:n_pathways, function(i) {
        size <- sample(20:100, 1)
        paste0("gene", sample(n_genes, size))
    })
    names(pathways) <- paste0("pathway", 1:n_pathways)

    # Should complete in reasonable time
    time <- system.time({
        result <- superfastGSEA(pathways, stats, minSize = 15)
    })

    # Should complete in under 30 seconds (conservative)
    expect_true(time["elapsed"] < 30)

    # Should return results for all valid pathways
    expect_true(nrow(result) > 0)
})

test_that("handles large gene sets", {
    skip_on_cran()

    set.seed(42)
    n_genes <- 10000

    stats <- setNames(rnorm(n_genes), paste0("gene", 1:n_genes))

    pathways <- list(
        large1 = paste0("gene", sample(n_genes, 500)),
        large2 = paste0("gene", sample(n_genes, 300))
    )

    # Should complete without timeout
    time <- system.time({
        result <- superfastGSEA(pathways, stats, minSize = 100, maxSize = 600)
    })

    expect_true(time["elapsed"] < 10)
    expect_equal(nrow(result), 2)
})
