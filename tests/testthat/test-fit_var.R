test_that("fit_var returns correct dimensions", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))
  result <- timecop:::fit_var(obj@gamma_hat, obj@Gamma_hat, obj@d)

  expect_equal(dim(result), c(2, 2))
})

test_that("fit_var with identity Gamma returns gamma", {
  d <- 2
  gamma_hat <- matrix(c(0.3, 0.1, 0.1, 0.3), d, d)
  Gamma_hat <- diag(d)

  result <- timecop:::fit_var(gamma_hat, Gamma_hat, d)

  expect_equal(result, gamma_hat, tolerance = 1e-10)
})
