# ---- constructor -----------------------------------------------------------

test_that("timecop_multitask constructs from a list of subjects", {
  d <- make_multitask_data()
  obj <- timecop_multitask(d$data, d$family)

  expect_s4_class(obj, "timecop_multitask")
  expect_equal(obj@K, d$K)
  expect_equal(obj@d, d$d)
  expect_length(obj@subjects, d$K)
  expect_true(all(vapply(obj@subjects, is, logical(1), class2 = "timecop")))
})

test_that("timecop_multitask errors on a non-list and reports bad subjects", {
  d <- make_multitask_data()
  expect_error(timecop_multitask(d$data[[1]], d$family), "non-empty list")

  bad <- d$data
  bad[[1]][1, 1] <- NA
  expect_error(timecop_multitask(bad, d$family), "subject 1")
})

# ---- temporal solver -------------------------------------------------------

test_that("multitask_pgd recovers per-subject OLS at lambda = 0", {
  set.seed(3); d <- 2; K <- 2
  B <- matrix(c(0.5, 0.1, 0.0, 0.4), 2, 2)
  mk <- function(n) {
    x <- matrix(0, d, n + 1)
    for (t in 2:(n + 1)) x[, t] <- B %*% x[, t - 1] + rnorm(d, 0, 0.3)
    list(G = crossprod(t(x[, 1:n])), H = crossprod(t(x[, 2:(n + 1)]), t(x[, 1:n])))
  }
  s <- lapply(seq_len(K), function(k) mk(2000))
  G <- lapply(s, `[[`, "G"); H <- lapply(s, function(z) t(z$H))
  f <- timecop:::multitask_pgd(G, H, 0, 0, max_iter = 3000, tol = 1e-10)
  ols <- lapply(seq_len(K), function(k) t(solve(G[[k]], H[[k]])))
  expect_lt(max(vapply(seq_len(K),
                       function(k) max(abs(f$B[[k]] - ols[[k]])), numeric(1))), 1e-4)
})

test_that("multitask_pgd with Omega = identity matches Omega = NULL", {
  set.seed(4); d <- 2; K <- 2
  G <- lapply(seq_len(K), function(k) { A <- matrix(rnorm(d * d), d, d); crossprod(A) + diag(d) })
  H <- lapply(seq_len(K), function(k) matrix(rnorm(d * d), d, d))
  f0 <- timecop:::multitask_pgd(G, H, 0.2, 0.2, max_iter = 1000, tol = 1e-10)
  I  <- lapply(seq_len(K), function(k) diag(d))
  f1 <- timecop:::multitask_pgd(G, H, 0.2, 0.2, Omega = I, max_iter = 1000, tol = 1e-10)
  expect_equal(f0$mu, f1$mu)
  expect_equal(f0$delta, f1$delta)
})

# ---- precision solver ------------------------------------------------------

test_that("multitask_glasso_admm returns PD precisions with Omega = M + E", {
  set.seed(5); d <- 3; K <- 2
  Th <- diag(1.5, d); Th[1, 2] <- Th[2, 1] <- -0.35; Th[2, 3] <- Th[3, 2] <- -0.30
  rmvn <- function(n, Sigma) matrix(rnorm(n * ncol(Sigma)), n) %*% chol(Sigma)
  S1 <- cov(rmvn(1500, solve(Th)))
  S2 <- cov(rmvn(1500, solve(Th)))
  fm <- timecop:::multitask_glasso_admm(list(S1, S2), 0.08, 0.08,
                                        weights = c(750, 750))

  expect_true(fm$converged)
  expect_true(all(vapply(fm$Omega,
    function(O) min(eigen(O, symmetric = TRUE, only.values = TRUE)$values) > 0,
    logical(1))))
  expect_true(all(abs(diag(fm$M)) < 1e-10))              # diag(M) == 0
  for (k in seq_len(K)) {
    expect_equal(fm$Omega[[k]], fm$M + fm$E[[k]], tolerance = 1e-3)
  }
})

# ---- ordinary multitask fit ------------------------------------------------

test_that("fit_multitask (gvar = FALSE) returns a multitask fit", {
  d <- make_multitask_data()
  obj <- timecop_multitask(d$data, d$family)
  fit <- fit_multitask(obj, n_lambda_mu = 4, n_lambda_delta = 4)

  expect_s3_class(fit, "timecop_multitask_fit")
  expect_true(is.matrix(fit$mu_hat))
  expect_length(fit$delta_hat, d$K)
  expect_true(all(is.finite(fit$mu_hat)))
  expect_true(all(is.finite(fit$ebic_grid)))
})

# ---- multitask graphical VAR (gvar = TRUE) ---------------------------------

test_that("fit_multitask (gvar = TRUE) fits the graphical VAR over the grid", {
  d <- make_multitask_data()
  obj <- timecop_multitask(d$data, d$family)
  fit <- fit_multitask(obj, gvar = TRUE,
                       n_lambda_mu = 3, n_lambda_delta = 3, n_lambda_prec = 3)

  expect_s3_class(fit, "timecop_multitask_gvar_fit")
  expect_length(fit$Omega_hat, d$K)
  # precision networks are positive definite
  expect_true(all(vapply(fit$Omega_hat,
    function(O) min(eigen(O, symmetric = TRUE, only.values = TRUE)$values) > 0,
    logical(1))))
  # diag(M) == 0 by construction
  expect_true(all(abs(diag(fit$M_hat)) < 1e-8))
  # selection returns the grid minimum
  expect_equal(dim(fit$ebic_grid), c(3, 3, 3, 3))
  expect_equal(fit$ebic, min(fit$ebic_grid))
})

test_that("gvar = TRUE with supplied scalar penalties collapses the grid", {
  d <- make_multitask_data()
  obj <- timecop_multitask(d$data, d$family)
  fit <- fit_multitask(obj, gvar = TRUE,
                       lambda_mu = 0.05, lambda_delta = 0.05,
                       lambda_M = 0.1, lambda_E = 0.1)

  expect_equal(dim(fit$ebic_grid), c(1L, 1L, 1L, 1L))
  expect_equal(fit$lambda_M, 0.1)
})

test_that("coordinate search matches the full grid selection", {
  d <- make_multitask_data()
  obj <- timecop_multitask(d$data, d$family)
  fg <- fit_multitask(obj, gvar = TRUE, search = "grid",
                      n_lambda_mu = 3, n_lambda_delta = 3, n_lambda_prec = 3)
  fc <- fit_multitask(obj, gvar = TRUE, search = "coordinate",
                      n_lambda_mu = 3, n_lambda_delta = 3, n_lambda_prec = 3)

  # same selected penalties and EBIC (coordinate visits a subset of the grid)
  expect_equal(fc$ebic, fg$ebic, tolerance = 1e-6)
  expect_equal(c(fc$lambda_mu, fc$lambda_delta, fc$lambda_M, fc$lambda_E),
               c(fg$lambda_mu, fg$lambda_delta, fg$lambda_M, fg$lambda_E))
  # coordinate evaluates fewer points (unvisited = NA), grid evaluates all
  expect_true(sum(is.na(fc$ebic_grid)) > 0)
  expect_true(all(is.finite(fg$ebic_grid)))
})

# ---- KKT anchoring of the automatic penalty grids --------------------------
# The top of each auto grid is the smallest penalty at which the corresponding
# block is fully sparse (subgradient condition at zero). Fitting AT the grid
# top must therefore return an empty network on every axis. This guards the
# calibration of all four grids (see multitask_precision_grid_diagnosis.md).

test_that("auto grid tops are KKT zeroing thresholds (fully sparse fits)", {
  d <- make_multitask_data()
  obj <- timecop_multitask(d$data, d$family)
  p <- obj@p; K <- obj@K; N <- obj@N
  S0  <- lapply(obj@subjects, function(s) s@cov_z_hat[,, p + 1L])
  S10 <- lapply(obj@subjects, function(s) s@cov_z_hat[,, p])
  G   <- lapply(seq_len(K), function(k) N[k] * S0[[k]])
  H   <- lapply(seq_len(K), function(k) N[k] * t(S10[[k]]))

  grids <- timecop:::multitask_lambda_grids(H, S0, S10, N, 3L, 3L, 3L,
                                            gvar = TRUE)

  # temporal: at (lam_mu_max, lam_delta_max) the temporal networks are empty
  tf <- timecop:::multitask_pgd(G, H,
                                grids$lambda_mu_seq[1], grids$lambda_delta_seq[1])
  expect_lt(sum(abs(tf$mu)), 1e-8)
  expect_lt(sum(vapply(tf$delta, function(D) sum(abs(D)), numeric(1))), 1e-8)

  # precision: at (lam_M_max, lam_E_max) the precision networks are empty
  B_yw  <- lapply(seq_len(K), function(k) S10[[k]] %*% solve(S0[[k]]))
  S_eps <- lapply(seq_len(K), function(k)
    timecop:::gvar_resid_cov(B_yw[[k]], S0[[k]], S0[[k]], t(S10[[k]]), S10[[k]]))
  f <- timecop:::multitask_glasso_admm(S_eps,
                                       grids$lambda_M_seq[1], grids$lambda_E_seq[1],
                                       weights = N / 2)
  expect_equal(sum(abs(f$M[upper.tri(f$M)]) > 1e-6), 0)
  expect_equal(sum(vapply(f$E, function(E) sum(abs(E[upper.tri(E)]) > 1e-6),
                          numeric(1))), 0)

  # and slightly below the anchors the penalty no longer fully sparsifies
  f2 <- timecop:::multitask_glasso_admm(S_eps,
                                        grids$lambda_M_seq[1] * 0.2,
                                        grids$lambda_E_seq[1] * 0.2,
                                        weights = N / 2)
  expect_gt(sum(abs(f2$M[upper.tri(f2$M)]) > 1e-6) +
              sum(vapply(f2$E, function(E) sum(abs(E[upper.tri(E)]) > 1e-6),
                         numeric(1))), 0)
})

# ---- PD guards --------------------------------------------------------------

test_that("pd_approx = TRUE on healthy data only standardizes the lag-0 diagonal", {
  d <- make_multitask_data()
  obj_F <- timecop_multitask(d$data, d$family, pd_approx = FALSE)
  obj_T <- timecop_multitask(d$data, d$family, pd_approx = TRUE)
  for (k in seq_len(d$K)) {
    p  <- obj_T@subjects[[k]]@p
    zT <- obj_T@subjects[[k]]@cov_z_hat
    zF <- obj_F@subjects[[k]]@cov_z_hat
    # lag-0 diagonal is set to exactly 1 (the true latent self-correlation)
    expect_identical(diag(zT[,, p + 1L]), rep(1, d$d))
    # everything else is untouched on the no-repair path
    off <- upper.tri(zT[,, p + 1L]) | lower.tri(zT[,, p + 1L])
    expect_identical(zT[,, p + 1L][off], zF[,, p + 1L][off])
    expect_identical(zT[,, p], zF[,, p])
  }
})

test_that("EBIC functions return Inf for non-PD residual covariances", {
  dd <- 3; K <- 2
  S0  <- replicate(K, diag(dd), simplify = FALSE)
  S10 <- replicate(K, diag(1.5, dd), simplify = FALSE)   # implies B = 1.5 I
  N   <- c(100, 100)
  mu  <- diag(1.5, dd)                                    # S_eps = -1.25 I: non-PD
  delta <- replicate(K, matrix(0, dd, dd), simplify = FALSE)
  B     <- replicate(K, diag(1.5, dd), simplify = FALSE)

  eb <- timecop:::multitask_ebic(mu, delta, B, S0, S10, N)
  expect_identical(eb$ebic, Inf)

  Omega <- replicate(K, diag(dd), simplify = FALSE)
  M     <- matrix(0, dd, dd)
  E     <- replicate(K, diag(dd), simplify = FALSE)
  eb2 <- timecop:::multitask_gvar_ebic(mu, delta, B, Omega, M, E, S0, S10, N)
  expect_identical(eb2$ebic, Inf)
})

test_that("gvar fit warns on non-PD joint latent covariance", {
  d <- make_multitask_data()
  obj <- timecop_multitask(d$data, d$family)
  # tamper one subject's lag-1 block so the joint matrix is clearly non-PSD
  s <- obj@subjects[[1]]
  s@cov_z_hat[,, s@p] <- matrix(0.99, d$d, d$d)
  obj@subjects[[1]] <- s
  # fires the upfront warning, and (with the single grid point invalid) the
  # all-Inf guard stops with actionable advice instead of returning garbage
  expect_warning(
    expect_error(
      fit_multitask(obj, gvar = TRUE,
                    lambda_mu = 10, lambda_delta = 10,
                    lambda_M = 1, lambda_E = 1, max_outer = 2L),
      "no valid model"
    ),
    "pd_approx = TRUE"
  )
})

test_that("adaptive-weighted auto grid tops are weighted KKT thresholds", {
  d <- make_multitask_data()
  obj <- timecop_multitask(d$data, d$family)
  p <- obj@p; K <- obj@K; N <- obj@N
  S0  <- lapply(obj@subjects, function(s) s@cov_z_hat[,, p + 1L])
  S10 <- lapply(obj@subjects, function(s) s@cov_z_hat[,, p])
  G   <- lapply(seq_len(K), function(k) N[k] * S0[[k]])
  H   <- lapply(seq_len(K), function(k) N[k] * t(S10[[k]]))

  aw <- timecop:::multitask_adaptive_weights(obj@subjects, S0, S10,
                                             source = "unpenalized")
  grids <- timecop:::multitask_lambda_grids(H, S0, S10, N, 3L, 3L, 3L,
                                            W_mu = aw$W_mu, W_delta = aw$W_delta,
                                            W_M = aw$W_M, W_E = aw$W_E,
                                            gvar = TRUE)

  # weighted anchors sit BELOW the uniform ones (all adaptive weights > 1 here)
  g0 <- timecop:::multitask_lambda_grids(H, S0, S10, N, 3L, 3L, 3L, gvar = TRUE)
  expect_lt(grids$lambda_M_seq[1], g0$lambda_M_seq[1])
  expect_lt(grids$lambda_E_seq[1], g0$lambda_E_seq[1])

  # temporal: at the weighted tops the weighted fit is fully sparse
  tf <- timecop:::multitask_pgd(G, H,
                                grids$lambda_mu_seq[1], grids$lambda_delta_seq[1],
                                W_mu = aw$W_mu, W_delta = aw$W_delta)
  expect_lt(sum(abs(tf$mu)), 1e-8)
  expect_lt(sum(vapply(tf$delta, function(D) sum(abs(D)), numeric(1))), 1e-8)

  # precision: at the weighted tops the weighted ADMM fit is fully sparse
  B_yw  <- lapply(seq_len(K), function(k) S10[[k]] %*% solve(S0[[k]]))
  S_eps <- lapply(seq_len(K), function(k)
    timecop:::gvar_resid_cov(B_yw[[k]], S0[[k]], S0[[k]], t(S10[[k]]), S10[[k]]))
  f <- timecop:::multitask_glasso_admm(S_eps,
                                       grids$lambda_M_seq[1], grids$lambda_E_seq[1],
                                       weights = N / 2,
                                       W_M = aw$W_M, W_E = aw$W_E)
  expect_equal(sum(abs(f$M[upper.tri(f$M)]) > 1e-6), 0)
  expect_equal(sum(vapply(f$E, function(E) sum(abs(E[upper.tri(E)]) > 1e-6),
                          numeric(1))), 0)
})
