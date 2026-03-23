# superfastGSEA

Super Fast Gene Set Enrichment Analysis with a Rust backend.

## Overview

`superfastGSEA` is an R package that provides high-performance Gene Set Enrichment Analysis (GSEA) using a Rust backend with Gamma distribution approximation. It achieves a 10–50x speedup compared to `fgsea` while maintaining statistical accuracy through anchor sampling and parallel computation via `rayon`.

## Installation

```r
# Install from source (requires Rust/Cargo)
devtools::install_github("your-org/superfastgsea")
```

## Requirements

- R >= 4.0
- Rust toolchain (cargo, rustc) — install via [rustup.rs](https://rustup.rs/)

## Usage

```r
library(superfastgsea)

# Run GSEA with DESeq2 statistics
results <- superfastGSEA(
  stats    = deseq2_stats,   # Named numeric vector of test statistics
  pathways = pathway_list    # Named list of gene sets
)

# The function will receive the stats vector and dispatch computation
# to the Rust core for parallel pathway scoring.
print(results)
```

## How It Works

1. **R front-end** — validates inputs, handles NAs, and passes data pointers (not copies) to Rust via `rextendr`.
2. **Rust core** — uses `rayon` for fearless parallelism across pathways; each worker computes ES, fits a Gamma distribution, and returns NES + p-value.
3. **Fallback** — pathways too small for Gamma fitting automatically fall back to exact permutation.

## Performance

| Package       | Time (1 000 pathways) | Notes                  |
|---------------|-----------------------|------------------------|
| fgsea         | ~10 s                 | R + C++ baseline       |
| superfastGSEA | ~0.3 s                | Rust + rayon parallel  |

## License

MIT
