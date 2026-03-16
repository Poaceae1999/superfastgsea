# superfastgsea

Super Fast Gene Set Enrichment Analysis with Rust Backend.

## Overview

`superfastgsea` is a high-performance R package for Gene Set Enrichment Analysis (GSEA) that uses a Rust backend with Gamma distribution approximation. It achieves 10-50x speedup compared to fgsea while maintaining statistical accuracy through anchor sampling and parallel computation via rayon.

The API is compatible with `fgsea::fgsea()`, so you can receive results in the same data.frame format and integrate seamlessly into existing workflows.

## Installation

```r
# Install from source (requires Rust toolchain)
# install.packages("devtools")
devtools::install_github("user/superfastgsea")
```

### System Requirements

- R (>= 4.0)
- Rust toolchain (`cargo` and `rustc`)
- `rextendr` (>= 0.4.2)

## Usage

```r
library(superfastgsea)

# Create example data
set.seed(42)
stats <- setNames(rnorm(1000), paste0("gene", 1:1000))
stats <- sort(stats, decreasing = TRUE)

pathways <- list(
  pathway1 = paste0("gene", 1:50),
  pathway2 = paste0("gene", 500:550)
)

# Fast analysis with Gamma method (default)
result <- superfastGSEA(pathways, stats, method = "gamma")

# More accurate analysis with Empirical+GPD method
result_gpd <- superfastGSEA(pathways, stats, method = "empirical_gpd")
```

## Methods

| Method | Description |
|--------|-------------|
| `gamma` | Gamma distribution approximation (fast, default) |
| `gamma_adaptive` | Gamma with adaptive sampling (stops when converged) |
| `gamma_two_pass` | Two-pass coarse-to-fine (fastest for large datasets) |
| `empirical_gpd` | Empirical CDF + GPD for tails (more accurate, uses EVT) |
| `hybrid` | Use fgsea for small pathways, Gamma for large (highest accuracy, requires fgsea) |

## Output

You receive a data.frame with the following columns:

| Column | Description |
|--------|-------------|
| `pathway` | Pathway name |
| `pval` | Nominal p-value |
| `padj` | Adjusted p-value |
| `ES` | Enrichment Score |
| `NES` | Normalized Enrichment Score |
| `size` | Number of genes in pathway after filtering |
| `leadingEdge` | List of leading edge gene names |

## License

MIT
