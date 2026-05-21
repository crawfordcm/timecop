# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Rules

- **Always ask the user for confirmation before editing any R files (R/*.R, tests/testthat/*.R).** Describe the planned changes and wait for approval before making edits.
- **Never add `Co-Authored-By: Claude` to commit messages.**

## Package Overview

`timecop` is an R package for estimating latent Gaussian VAR (Vector Autoregressive) models for discrete-valued multivariate time series. It supports Bernoulli, Poisson, and Gaussian marginal distributions through a latent variable approach with Hermite polynomial expansion. Currently only VAR order p=1 is supported.

## Build & Development Commands

This is a standard R package using devtools. From R:

```r
devtools::load_all()       # Load package for interactive development
devtools::document()       # Regenerate roxygen2 documentation (man/ and NAMESPACE)
devtools::check()          # Full R CMD check
devtools::build()          # Build package tarball
devtools::install()        # Install locally
```

From the command line:
```bash
R CMD build .
R CMD check timecop_0.0.0.9000.tar.gz
R CMD INSTALL .
```

Documentation uses roxygen2 with markdown enabled. The man/ directory is auto-generated — edit roxygen comments in R/ files, not man/*.Rd files directly.

There is no test suite (no testthat or tests/ directory).

## Architecture

### S4 Class System

The package uses R's S4 class system (via `methods`). The central class is `timecop` defined in `R/timecopObjectClass.R` with these key slots:

- `data` — d×n matrix (variables × time points; user provides n×d, constructor transposes)
- `cov_x_hat` — observed covariance arrays at multiple lags
- `cov_z_hat` — latent Gaussian covariance arrays (estimated via link functions)
- `gamma_hat` / `Gamma_hat` — stacked lag covariances and Toeplitz matrix for Yule-Walker estimation
- `family` — list of marginal distributions per variable

### Public API (3 exported functions)

1. **`timecop()`** — Constructor. Takes n×d data matrix and family list, computes observed covariances, applies link/inverse-link functions to get latent covariances, builds Yule-Walker matrices.
2. **`fit_timecop()`** — S4 generic method. Runs Yule-Walker estimation (`fit_var`) and computes standard errors (`se_var`). Returns a list of `(estimates, standard_errors)`.
3. **`latent_var_sim()`** — Simulates multivariate time series from a latent Gaussian VAR with discrete marginals.

### Estimation Pipeline

The constructor `timecop()` runs this pipeline:

1. `setup_data()` → validate and transpose input to d×n
2. `check_marginals()` → validate family distributions
3. `observed_var_cov()` → compute observed covariance/correlation matrices at lags 0..2p
4. `latent_var_link()` → compute link function coefficients via `link_coefs()` → `hermite_coefs()` (uses pre-computed Hermite polynomials from `sysdata.rda`)
5. `latent_var_invlink()` → invert link via `interpolation()` → `nat_spline()` (natural cubic spline)
6. Build `gamma_hat` (stacked lag covariances) and `Gamma_hat` (Toeplitz matrix)
7. Optionally `check_pd()` → ensure positive definiteness

Then `fit_timecop()` runs:

1. `fit_var()` → Yule-Walker: `solve(Gamma_hat) %*% t(gamma_hat)`
2. `se_var()` → standard errors via `numderiv()` (Jacobian) and `longrun_var()` (long-run variance)

### Matrix Algebra Utilities

`R/matrixcalc.R` implements Magnus & Neudecker matrix differentiation operators: `vec()`, `vech()`, `vecp()`, `vecd()`, commutation matrices, duplication matrices, elimination matrices. These are used for computing numerical derivatives in the standard error calculation.

### Pre-computed Data

`R/sysdata.rda` contains a `Polys` object with pre-computed Hermite polynomials used by `hermite_coefs()`.

## Key Dependencies

- `mvtnorm` — multivariate normal distributions
- `numDeriv` — numerical derivatives (Jacobian computation)
- `polynom` — polynomial operations for Hermite coefficients
- `cointReg` — long-run variance estimation (`getLongRunVar`)
- `Matrix` — `nearPD` for positive-definite matrix approximation
- `expm` — matrix exponentials (used in simulation)
