test_that("compact ordinal levels expand to align with all variables", {
  data <- rbind(
    c(0, 1, 0, 1, 0, 1),
    c(1, 2, 3, 1, 2, 3),
    c(-1, 0, 1, -1, 0, 1),
    c(5, 6, 7, 8, 5, 6)
  )
  family <- list("Bernoulli", "Ordinal", "Gaussian", "Ordinal")

  result <- timecop:::setup_ordinal_levels(
    data = data,
    family = family,
    ordinal_levels = list(1:3, 5:8),
    d = 4
  )

  expect_equal(result, list(NULL, 1:3, NULL, 5:8))
})

test_that("invalid ordinal levels are rejected", {
  expect_error(
    timecop:::setup_ordinal_levels(
      data = matrix(c(0, 1, 0, 1), nrow = 1),
      family = list("Ordinal"),
      ordinal_levels = list(0:1),
      d = 1
    ),
    "Bernoulli"
  )

  expect_error(
    timecop:::setup_ordinal_levels(
      data = matrix(c(1, 2, 3, 1), nrow = 1),
      family = list("Ordinal"),
      ordinal_levels = list(c(1, 2, 5)),
      d = 1
    ),
    "consecutive"
  )

  expect_error(
    timecop:::setup_ordinal_levels(
      data = matrix(c(1, 2, 3, 1), nrow = 1),
      family = list("Ordinal"),
      ordinal_levels = list(c(1, 3, 2)),
      d = 1
    ),
    "strictly increasing"
  )

  expect_error(
    timecop:::setup_ordinal_levels(
      data = matrix(c(1, 2, 3, 4), nrow = 1),
      family = list("Ordinal"),
      ordinal_levels = list(1:3),
      d = 1
    ),
    "not listed"
  )
})

test_that("ordinal marginal probabilities follow supplied level order", {
  data <- matrix(c(1, 1, 2, 2, 2, 3, 3, 3), nrow = 1)

  result <- timecop:::estimate_marginal_params(
    data = data,
    family = list("Ordinal"),
    ordinal_levels = list(1:3),
    d = 1
  )

  expect_equal(result[[1]], c(2 / 8, 3 / 8, 3 / 8))

  expect_error(
    timecop:::estimate_marginal_params(
      data = data,
      family = list("Ordinal"),
      ordinal_levels = list(1:4),
      d = 1
    ),
    "not observed"
  )
})

test_that("three- and four-category ordinal variables can be estimated", {
  set.seed(20261002)

  sim <- latent_var_sim(
    d = 2,
    n = 500,
    p = 1,
    param = list(
      c(0.20, 0.50, 0.30),
      c(0.10, 0.20, 0.30, 0.40)
    ),
    phi_lv = matrix(
      c(0.35, 0.10,
        0.15, 0.30),
      nrow = 2,
      byrow = TRUE
    ),
    family = list("Ordinal", "Ordinal")
  )

  expect_equal(sim$ordinal_levels, list(0:2, 0:3))
  expect_true(all(sim$X_t[1, ] %in% 0:2))
  expect_true(all(sim$X_t[2, ] %in% 0:3))

  shifted_data <- sim$X_t
  shifted_data[1, ] <- shifted_data[1, ] + 1
  shifted_data[2, ] <- shifted_data[2, ] + 5

  obj <- timecop(
    data = t(shifted_data),
    family = list("Ordinal", "Ordinal"),
    ordinal_levels = list(1:3, 5:8)
  )

  estimates <- timecop:::fit_var(
    gamma_hat = obj@gamma_hat,
    Gamma_hat = obj@Gamma_hat,
    d = obj@d
  )

  expect_equal(obj@ordinal_levels, list(1:3, 5:8))
  expect_equal(obj@marg_num, 5)
  expect_equal(dim(estimates), c(2, 2))
  expect_true(all(is.finite(estimates)))
})

test_that("two-category ordinal simulation recommends Bernoulli", {
  expect_error(
    latent_var_sim(
      d = 1,
      n = 100,
      p = 1,
      param = list(c(0.5, 0.5)),
      phi_lv = matrix(0.3, 1, 1),
      family = list("Ordinal")
    ),
    "Bernoulli"
  )
})
