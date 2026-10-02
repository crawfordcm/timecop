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

test_that("interpolation handles values on the boundary knots without error", {
  coef <- c(0.8, 0, 0.1, rep(0, 97))
  u <- seq(-1, 1, length.out = 21)
  pow <- seq_along(coef)
  knot <- vapply(u, function(uu) sum(coef * uu^pow), numeric(1))
  n <- length(knot)

  # a value exactly on the first knot must invert to u[1], not error
  expect_equal(timecop:::interpolation(coef, u, knot[1]), u[1])
  # and on the last knot to u[n]
  expect_equal(timecop:::interpolation(coef, u, knot[n]), u[n])
  # interior knots are exact too
  expect_equal(timecop:::interpolation(coef, u, knot[5]), u[5])
})
