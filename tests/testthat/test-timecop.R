test_that("constructor creates valid S4 object with Bernoulli", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))
  expect_s4_class(obj, "timecop")
  expect_equal(obj@d, 2)
  expect_equal(obj@n, 200)
  expect_equal(obj@p, 1)
})

test_that("constructor works with Poisson", {
  sim <- make_test_sim_poisson()
  obj <- timecop(data = t(sim$X_t), family = list("Poisson", "Poisson"))
  expect_s4_class(obj, "timecop")
})

test_that("constructor rejects NULL data", {
  expect_error(timecop(data = NULL, family = list("Bernoulli")),
               "not supplied")
})

test_that("constructor rejects NULL family", {
  sim <- make_test_sim()
  expect_error(timecop(data = t(sim$X_t), family = NULL),
               "marginal distributions")
})

test_that("constructor rejects unsupported family", {
  sim <- make_test_sim()
  expect_error(timecop(data = t(sim$X_t), family = list("NegBin", "Bernoulli")),
               "not supported")
})

test_that("constructor rejects p != 1", {
  sim <- make_test_sim()
  expect_error(timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"), p = 2),
               "Only lag 1")
})

test_that("constructor rejects data with NAs", {
  sim <- make_test_sim()
  bad_data <- t(sim$X_t)
  bad_data[1, 1] <- NA
  expect_error(timecop(data = bad_data, family = list("Bernoulli", "Bernoulli")),
               "NA values")
})

test_that("constructor works with pd_approx = TRUE", {
  sim <- make_test_sim()
  obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"),
                 pd_approx = TRUE)
  expect_s4_class(obj, "timecop")
  expect_true(obj@pd_approx)
})
