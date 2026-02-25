# timecop

Estimation of latent Gaussian VAR(1) models for discrete-valued multivariate
time series. Supports Bernoulli, Poisson, and Gaussian marginal distributions
through a latent variable approach using Hermite polynomial expansion.

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("crawfordcm/timecop")
```

## Usage

```r
library(timecop)

# Simulate bivariate Bernoulli time series
sim <- latent_var_sim(
  d = 2, n = 200, p = 1,
  param = list(0.5, 0.5),
  phi_lv = matrix(c(0.4, 0.2, 0.2, 0.4), 2, 2),
  family = list("Bernoulli", "Bernoulli")
)

# Construct timecop object
obj <- timecop(
  data = t(sim$X_t),
  family = list("Bernoulli", "Bernoulli")
)

# Fit VAR(1) model
results <- fit_timecop(obj)
results$estimates       # coefficient matrix
results$se              # standard errors
summary(results)        # formatted summary table
```

## License

MIT
