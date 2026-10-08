test_that("se_var returns correct dimensions", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))

  se <- timecop:::se_var(
    data = obj@data,
    gamma_hat = obj@gamma_hat,
    Gamma_hat = obj@Gamma_hat,
    cov_x_hat = obj@cov_x_hat,
    d = obj@d,
    p = obj@p,
    n = obj@n,
    family = obj@family,
    ordinal_levels = obj@ordinal_levels,
    marg_num = obj@marg_num,
    corr = obj@corr
  )

  expect_equal(dim(se), c(2, 2))
})

test_that("se_var returns positive finite values", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))

  se <- timecop:::se_var(
    data = obj@data,
    gamma_hat = obj@gamma_hat,
    Gamma_hat = obj@Gamma_hat,
    cov_x_hat = obj@cov_x_hat,
    d = obj@d,
    p = obj@p,
    n = obj@n,
    family = obj@family,
    ordinal_levels = obj@ordinal_levels,
    marg_num = obj@marg_num,
    corr = obj@corr
  )

  expect_true(all(se > 0))
  expect_true(all(is.finite(se)))
})
