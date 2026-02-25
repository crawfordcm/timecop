test_that("se_var returns correct dimensions", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))

  se <- timecop:::se_var(
    obj@data, obj@gamma_hat, obj@Gamma_hat, obj@cov_x_hat,
    obj@d, obj@p, obj@n, obj@family, obj@marg_num, obj@corr
  )

  expect_equal(dim(se), c(2, 2))
})

test_that("se_var returns positive finite values", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))

  se <- timecop:::se_var(
    obj@data, obj@gamma_hat, obj@Gamma_hat, obj@cov_x_hat,
    obj@d, obj@p, obj@n, obj@family, obj@marg_num, obj@corr
  )

  expect_true(all(se > 0))
  expect_true(all(is.finite(se)))
})
