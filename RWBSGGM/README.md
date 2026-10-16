# RWBSGGM

R implementation accompanying the manuscript

**Efficient Oracle-based Estimation and Inference for Multiple Change Points in High-Dimensional Gaussian Graphical Models**

by Xianru Wang, Zhi Yang, Bin Liu, and Yufeng Liu.

## Package overview

This package implements the RWBS procedure used in the paper, including:

- initial multiple change-point estimation by refined wild binary segmentation;
- dynamic-programming refinement;
- graphical-model estimation;
- confidence interval construction for change-point locations;
- auxiliary graph-generation and evaluation functions.

The package is organized so that the core functions are loaded automatically after installation. The simulation scripts are provided separately under `inst/scripts/`.

## Installation

The archive contains a standard R package source directory. After unzipping the
archive, install the package from a terminal with:

```bash
R CMD INSTALL RWBSGGM
```

Before installation, if the required packages are not already available, run in R:

```r
install.packages(c("glmnet", "MASS", "foreach", "doParallel"))
```

Alternatively, in R/RStudio, the unpacked source directory can be installed with
a local-package installer such as `remotes::install_local()`.

The core package imports `glmnet`; the simulation scripts additionally use
`MASS`, `foreach`, and `doParallel`.

## Quick example

After installation:

```r
library(RWBSGGM)

set.seed(123)
p <- 10
n <- 80

G <- chain.graph2(p, alpha = 2)
X1 <- MASS::mvrnorm(n/2, rep(0, p), solve(G$graph1))
X2 <- MASS::mvrnorm(n/2, rep(0, p), solve(G$graph2))
X <- rbind(X1, X2)

fit <- rwbs_fit(
  x = X,
  M = 5,
  delta = ceiling(0.1 * n),
  delta1 = 0.1,
  kmax = 5,
  do_inference = FALSE
)

fit$alpha.initial
fit$alpha.refine
```

The same example is available as:

```r
source(system.file("examples", "example_basic.R", package = "RWBSGGM"))
```

## Inference

To construct the confidence intervals used in the paper, supply the package's
Brownian-argmax quantile approximation:

```r
fit <- rwbs_fit(
  x = X,
  M = 30,
  delta = ceiling(0.1 * n),
  delta1 = 0.1,
  kmax = 10,
  do_inference = TRUE,
  argmax_quantile_fun = argmax.quantile
)

fit$confidence
```

## Reproducing simulation studies

From the unpacked source directory, the main simulation can be run by:

```bash
Rscript inst/scripts/run_RWBS.R 600 20 4
```

The arguments are:

1. `n`: sample size;
2. `p`: data dimension;
3. `alpha0`: number of nodes/strength parameter used in the graph-generating configuration.

The script performs the repeated simulation and writes CSV summaries to the current working directory.

Runtime experiments can be run by:

```bash
Rscript inst/scripts/run_time_rwbs.R 600 20
```

The scripts use the same core functions as the installed package and do not require manually sourcing files from `R/`.

## File structure

```text
RWBSGGM/
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── R/
│   ├── rwbs_core_functions.R
│   ├── graph_model_generators.R
│   ├── argmax.quantile.R
│   ├── helpers.nodewise.R
│   ├── has-standard.R
│   ├── standardize.R
│   └── ...
└── inst/
    ├── examples/
    │   └── example_basic.R
    └── scripts/
        ├── run_RWBS.R
        └── run_time_rwbs.R
```

## Reproducibility notes

The simulation scripts set a separate random seed for each replication. Parallel
workers are used for the repeated simulation study. Output files contain the
estimated change-point locations, estimated numbers of change points, Hausdorff
distances, confidence intervals, and related quantities used to summarize the
numerical results in the paper.
