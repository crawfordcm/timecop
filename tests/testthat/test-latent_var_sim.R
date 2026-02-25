test_that("latent_var_sim returns correct dimensions", {
  sim <- make_test_sim(n = 100)
  expect_equal(dim(sim$X_t), c(2, 100))
  expect_equal(dim(sim$Z_t), c(2, 100))
  expect_equal(dim(sim$A_true), c(2, 2, 1))
  expect_equal(dim(sim$Sigma_true), c(2, 2))
})

test_that("Bernoulli produces 0/1 values", {
  sim <- make_test_sim(n = 100)
  expect_true(all(sim$X_t %in% c(0, 1)))
})

test_that("Poisson produces non-negative integers", {
  sim <- make_test_sim_poisson(n = 100)
  expect_true(all(sim$X_t >= 0))
  expect_true(all(sim$X_t == floor(sim$X_t)))
})

test_that("latent_var_sim is reproducible with same seed", {
  sim1 <- make_test_sim(n = 50)
  sim2 <- make_test_sim(n = 50)
  expect_identical(sim1$X_t, sim2$X_t)
})

test_that("univariate case works", {
  set.seed(42)
  sim <- latent_var_sim(
    d = 1, n = 100, p = 1,
    param = list(0.5),
    phi_lv = matrix(0.3, 1, 1),
    family = list("Bernoulli")
  )
  expect_equal(dim(sim$X_t), c(1, 100))
  expect_true(all(sim$X_t %in% c(0, 1)))
})

test_that("latent_var_sim rejects invalid inputs", {
  expect_error(latent_var_sim(d = 0, n = 100, p = 1,
    param = list(0.5), phi_lv = matrix(0.3), family = list("Bernoulli")),
    "positive integer")
  expect_error(latent_var_sim(d = 2, n = 100, p = 2,
    param = list(0.5, 0.5), phi_lv = matrix(0, 2, 2),
    family = list("Bernoulli", "Bernoulli")),
    "p=1")
  expect_error(latent_var_sim(d = 2, n = 100, p = 1,
    param = list(0.5), phi_lv = matrix(0, 2, 2),
    family = list("Bernoulli", "Bernoulli")),
    "length d")
})
