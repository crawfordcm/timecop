# Golden-file regression tests: lock in point estimates and standard errors.
# Any change to estimation code that shifts values beyond 1e-6 will fail here.

test_that("Bernoulli estimates and SEs are unchanged (timecop_example)", {
  obj <- timecop(data = timecop_example, family = list("Bernoulli", "Bernoulli"))
  fit <- fit_timecop(obj)

  expected_est <- matrix(
    c(0.1890927, 0.2805873, 0.2816084, 0.4334893),
    nrow = 2, ncol = 2
  )
  expected_se <- matrix(
    c(0.1118608, 0.1078071, 0.11048718, 0.09476625),
    nrow = 2, ncol = 2
  )

  expect_equal(fit$estimates, expected_est, tolerance = 1e-6)
  expect_equal(fit$se,        expected_se,  tolerance = 1e-6)
})

test_that("Poisson estimates and SEs are unchanged (seed 42)", {
  set.seed(42)
  sim <- latent_var_sim(
    d = 2, n = 200, p = 1,
    param   = list(3, 5),
    phi_lv  = matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2),
    family  = list("Poisson", "Poisson")
  )
  obj <- timecop(data = t(sim$X_t), family = list("Poisson", "Poisson"))
  fit <- fit_timecop(obj)

  expected_est <- matrix(
    c(0.2154396, 0.1442337, 0.03062569, 0.34054859),
    nrow = 2, ncol = 2
  )
  expected_se <- matrix(
    c(0.09587402, 0.06273840, 0.06100206, 0.06101275),
    nrow = 2, ncol = 2
  )

  expect_equal(fit$estimates, expected_est, tolerance = 1e-6)
  expect_equal(fit$se,        expected_se,  tolerance = 1e-6)
})

test_that("Mixed Bernoulli+Gaussian estimates and SEs are unchanged (seed 7)", {
  set.seed(7)
  sim <- latent_var_sim(
    d = 2, n = 200, p = 1,
    param   = list(0.5, c(0, 1)),
    phi_lv  = matrix(c(0.4, 0.2, 0.2, 0.4), 2, 2),
    family  = list("Bernoulli", "Gaussian")
  )
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Gaussian"))
  fit <- fit_timecop(obj)

  expected_est <- matrix(
    c(0.5165779, 0.1992202, 0.09385076, 0.36556702),
    nrow = 2, ncol = 2
  )
  expected_se <- matrix(
    c(0.09445984, 0.08984790, 0.08709313, 0.08277307),
    nrow = 2, ncol = 2
  )

  expect_equal(fit$estimates, expected_est, tolerance = 1e-6)
  expect_equal(fit$se,        expected_se,  tolerance = 1e-6)
})

test_that("Mixed Bernoulli+Poisson estimates and SEs are unchanged (seed 99)", {
  set.seed(99)
  sim <- latent_var_sim(
    d = 2, n = 200, p = 1,
    param   = list(0.5, 3),
    phi_lv  = matrix(c(0.4, 0.2, 0.2, 0.4), 2, 2),
    family  = list("Bernoulli", "Poisson")
  )
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Poisson"))
  fit <- fit_timecop(obj)

  expected_est <- matrix(
    c(0.2599113, 0.1922599, 0.1640576, 0.3451767),
    nrow = 2, ncol = 2
  )
  expected_se <- matrix(
    c(0.10662377, 0.08966886, 0.08630834, 0.07981655),
    nrow = 2, ncol = 2
  )

  expect_equal(fit$estimates, expected_est, tolerance = 1e-6)
  expect_equal(fit$se,        expected_se,  tolerance = 1e-6)
})
