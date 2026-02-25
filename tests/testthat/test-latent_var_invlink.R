test_that("latent_var_invlink returns correct dimensions", {
  sim <- make_test_sim()
  data <- sim$X_t
  d <- nrow(data)
  n <- ncol(data)
  p <- 1
  k <- 100
  family <- list("Bernoulli", "Bernoulli")

  cov_x_hat <- timecop:::observed_var_cov(data, d, p, n, corr = FALSE)
  ell_ij_hat <- timecop:::latent_var_link(data, d, n, k, family, corr = FALSE)
  result <- timecop:::latent_var_invlink(cov_x_hat, d, p, ell_ij_hat)

  expect_equal(dim(result), dim(cov_x_hat))
})

test_that("latent_var_invlink values are finite", {
  sim <- make_test_sim()
  data <- sim$X_t
  d <- nrow(data)
  n <- ncol(data)
  p <- 1
  k <- 100
  family <- list("Bernoulli", "Bernoulli")

  cov_x_hat <- timecop:::observed_var_cov(data, d, p, n, corr = FALSE)
  ell_ij_hat <- timecop:::latent_var_link(data, d, n, k, family, corr = FALSE)
  result <- timecop:::latent_var_invlink(cov_x_hat, d, p, ell_ij_hat)

  expect_true(all(is.finite(result)))
})
