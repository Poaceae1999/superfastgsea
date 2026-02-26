# API Compatibility Tests
# Ensures superfastGSEA API matches fgsea

test_that("superfastGSEA accepts same input formats as fgsea", {
    # Create test data
    set.seed(42)
    stats <- setNames(rnorm(1000), paste0("gene", 1:1000))

    pathways <- list(
        pathway1 = paste0("gene", 1:50),
        pathway2 = paste0("gene", 500:550),
        pathway3 = paste0("gene", 950:1000)
    )

    # Should run without error
    expect_no_error(superfastGSEA(pathways, stats))
})

test_that("output format matches fgsea", {
    set.seed(42)
    stats <- setNames(rnorm(1000), paste0("gene", 1:1000))

    pathways <- list(
        pathway1 = paste0("gene", 1:50),
        pathway2 = paste0("gene", 500:550)
    )

    result <- superfastGSEA(pathways, stats)

    # Check required columns
    expected_cols <- c("pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge")
    expect_true(all(expected_cols %in% colnames(result)))

    # Check types
    expect_type(result$pathway, "character")
    expect_type(result$pval, "double")
    expect_type(result$padj, "double")
    expect_type(result$ES, "double")
    expect_type(result$NES, "double")
    expect_type(result$size, "integer")
    expect_type(result$leadingEdge, "list")

    # Check value ranges
    expect_true(all(result$pval >= 0 & result$pval <= 1))
    expect_true(all(result$padj >= 0 & result$padj <= 1))
    expect_true(all(result$ES >= -1 & result$ES <= 1))
})

test_that("minSize and maxSize filters work", {
    set.seed(42)
    stats <- setNames(rnorm(1000), paste0("gene", 1:1000))

    pathways <- list(
        small = paste0("gene", 1:5),      # 5 genes
        medium = paste0("gene", 10:50),   # 41 genes
        large = paste0("gene", 100:300)   # 201 genes
    )

    # Filter by size
    result <- superfastGSEA(pathways, stats, minSize = 10, maxSize = 100)

    # Only medium pathway should pass
    expect_equal(nrow(result), 1)
    expect_equal(result$pathway, "medium")
})

test_that("scoreType parameter works", {
    set.seed(42)
    # Create data where top genes are enriched
    stats <- setNames(c(rnorm(50, mean = 2), rnorm(950)), paste0("gene", 1:1000))
    stats <- sort(stats, decreasing = TRUE)

    pathways <- list(top_pathway = names(stats)[1:30])

    # Standard (two-tailed)
    result_std <- superfastGSEA(pathways, stats, scoreType = "std")
    expect_true(result_std$ES[1] > 0)

    # Positive only
    result_pos <- superfastGSEA(pathways, stats, scoreType = "pos")
    expect_true(result_pos$ES[1] >= 0)
})

test_that("handles NA values with warning", {
    stats <- setNames(c(1, 2, NA, 4, 5), paste0("gene", 1:5))
    pathways <- list(test = c("gene1", "gene2", "gene4"))

    expect_warning(
        superfastGSEA(pathways, stats),
        "Removed.*genes with NA"
    )
})

test_that("handles genes not in stats", {
    stats <- setNames(rnorm(100), paste0("gene", 1:100))
    pathways <- list(
        test = c("gene1", "gene2", "gene999", "geneXYZ")  # gene999 and geneXYZ don't exist
    )

    result <- superfastGSEA(pathways, stats, minSize = 1)
    # Should still work, just with reduced size
    expect_equal(result$size[1], 2)  # Only gene1 and gene2 found
})

test_that("empty pathways handled gracefully", {
    stats <- setNames(rnorm(100), paste0("gene", 1:100))
    pathways <- list()

    result <- superfastGSEA(pathways, stats)
    expect_equal(nrow(result), 0)
})

test_that("leadingEdge contains gene names", {
    set.seed(42)
    stats <- setNames(rnorm(100), paste0("gene", 1:100))
    stats <- sort(stats, decreasing = TRUE)

    pathways <- list(test = names(stats)[1:20])

    result <- superfastGSEA(pathways, stats, minSize = 1)

    # Leading edge should contain gene names, not indices
    expect_true(all(grepl("^gene", result$leadingEdge[[1]])))
})
