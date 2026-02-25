test_that("latent_var_link returns correct dimensions", {
  sim <- make_test_sim()
  data <- sim$X_t  # d x n
  d <- nrow(data)
  n <- ncol(data)
  k <- 100
  family <- list("Bernoulli", "Bernoulli")

  result <- timecop:::latent_var_link(data, d, n, k, family, corr = FALSE)

  expect_equal(dim(result), c(k, d, d))
})

test_that("latent_var_link values are finite", {
  sim <- make_test_sim()
  data <- sim$X_t
  d <- nrow(data)
  n <- ncol(data)
  k <- 100
  family <- list("Bernoulli", "Bernoulli")

  result <- timecop:::latent_var_link(data, d, n, k, family, corr = FALSE)

  expect_true(all(is.finite(result)))
})
