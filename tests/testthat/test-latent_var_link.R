test_that("latent_var_link returns correct dimensions", {
  sim <- make_test_sim()
  data <- sim$X_t  # d x n
  d <- nrow(data)
  n <- ncol(data)
  k <- 100
  family <- list("Bernoulli", "Bernoulli")
  ordinal_levels <- vector("list", d)

  result <- timecop:::latent_var_link(
    data = data,
    d = d,
    n = n,
    k = k,
    family = family,
    ordinal_levels = ordinal_levels,
    corr = FALSE
  )

  expect_equal(dim(result), c(k, d, d))
})

test_that("latent_var_link values are finite", {
  sim <- make_test_sim()
  data <- sim$X_t
  d <- nrow(data)
  n <- ncol(data)
  k <- 100
  family <- list("Bernoulli", "Bernoulli")
  ordinal_levels <- vector("list", d)

  result <- timecop:::latent_var_link(
    data = data,
    d = d,
    n = n,
    k = k,
    family = family,
    ordinal_levels = ordinal_levels,
    corr = FALSE
  )

  expect_true(all(is.finite(result)))
})
