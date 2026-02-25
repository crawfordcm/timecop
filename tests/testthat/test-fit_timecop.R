test_that("fit_timecop returns timecop_fit object", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))
  results <- fit_timecop(obj)

  expect_s3_class(results, "timecop_fit")
  expect_type(results, "list")
  expect_length(results, 3)
  expect_true(is.matrix(results$estimates))
  expect_true(is.matrix(results$se))
  expect_s4_class(results$obj, "timecop")
})

test_that("estimates are finite", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))
  results <- fit_timecop(obj)

  expect_true(all(is.finite(results$estimates)))
})

test_that("standard errors are positive", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))
  results <- fit_timecop(obj)

  expect_true(all(results$se > 0))
})

test_that("summary produces coefficient tables", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))
  results <- fit_timecop(obj)
  s <- summary(results)

  expect_s3_class(s, "summary.timecop_fit")
  expect_length(s$coefficients, 2)
  expect_equal(ncol(s$coefficients[[1]]), 4)
  expect_equal(nrow(s$coefficients[[1]]), 2)
  expect_true(all(s$coefficients[[1]][, "Pr(>|z|)"] >= 0))
  expect_true(all(s$coefficients[[1]][, "Pr(>|z|)"] <= 1))
})
