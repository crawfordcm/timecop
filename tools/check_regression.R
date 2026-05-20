#!/usr/bin/env Rscript
#
# check_regression.R — Numerical regression check for timecop
#
# Captures the package's numerical outputs across several representative
# scenarios and either saves them as a reference snapshot or compares the
# current outputs against a previously saved snapshot.
#
# Usage (from the package root):
#
#   Rscript tools/check_regression.R --generate
#       Run all scenarios, save outputs to tools/regression_reference.rds,
#       and print a summary of what was captured. Run this once on the current
#       codebase to establish the baseline.
#
#   Rscript tools/check_regression.R
#       Run all scenarios, compare against the saved reference, and report
#       any discrepancies. Exit code 0 = all pass, 1 = failures detected.
#
# The reference file should be committed to version control so that all
# collaborators share the same numerical baseline.
#
# Tolerance:
#   Comparisons use all.equal() with tolerance = 1e-6. Tighten or loosen
#   TOLERANCE below if needed for your platform.
#

REFERENCE_FILE <- file.path("tools", "regression_reference.rds")
TOLERANCE      <- 1e-6

# ── Helpers ────────────────────────────────────────────────────────────────────

cat_rule <- function(char = "-", width = 70) cat(strrep(char, width), "\n")

pass_fail <- function(ok, label) {
  status <- if (ok) "\033[32mPASS\033[0m" else "\033[31mFAIL\033[0m"
  cat(sprintf("  [%s] %s\n", status, label))
  ok
}

# Compare two numeric objects with all.equal(); return TRUE/FALSE + message.
num_equal <- function(current, reference, label, tol = TOLERANCE) {
  result <- all.equal(current, reference, tolerance = tol,
                      check.attributes = FALSE)
  ok <- isTRUE(result)
  if (!ok) {
    pass_fail(FALSE, label)
    for (msg in result) cat(sprintf("         %s\n", msg))
  } else {
    pass_fail(TRUE, label)
  }
  ok
}

# ── Package loading ─────────────────────────────────────────────────────────────

load_timecop <- function() {
  if (requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(devtools::load_all(".", quiet = TRUE))
    cat("Loaded timecop via devtools::load_all()\n")
  } else if (requireNamespace("timecop", quietly = TRUE)) {
    library(timecop)
    cat(sprintf("Loaded timecop %s from library\n", packageVersion("timecop")))
  } else {
    stop("Cannot load timecop. Install devtools or the package itself.")
  }
}

# ── Scenario definitions ────────────────────────────────────────────────────────
#
# Each scenario is a named list:
#   seed      — RNG seed for reproducibility
#   sim_args  — arguments to latent_var_sim()
#   tc_family — family argument to timecop() (must match sim_args$family)
#   label     — human-readable name shown in output

SCENARIOS <- list(

  bernoulli_2var = list(
    label     = "Bernoulli 2-variable",
    seed      = 42L,
    sim_args  = list(
      d       = 2,
      n       = 200,
      p       = 1,
      param   = list(0.5, 0.5),
      phi_lv  = matrix(c(0.4, 0.2, 0.2, 0.4), 2, 2),
      family  = list("Bernoulli", "Bernoulli")
    ),
    tc_family = list("Bernoulli", "Bernoulli")
  ),

  poisson_2var = list(
    label     = "Poisson 2-variable",
    seed      = 123L,
    sim_args  = list(
      d       = 2,
      n       = 200,
      p       = 1,
      param   = list(3, 5),
      phi_lv  = matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2),
      family  = list("Poisson", "Poisson")
    ),
    tc_family = list("Poisson", "Poisson")
  ),

  mixed_bern_pois = list(
    label     = "Mixed Bernoulli + Poisson",
    seed      = 7L,
    sim_args  = list(
      d       = 2,
      n       = 300,
      p       = 1,
      param   = list(0.4, 4),
      phi_lv  = matrix(c(0.35, 0.1, 0.15, 0.25), 2, 2),
      family  = list("Bernoulli", "Poisson")
    ),
    tc_family = list("Bernoulli", "Poisson")
  ),

  gaussian_2var = list(
    label     = "Gaussian 2-variable",
    seed      = 99L,
    sim_args  = list(
      d       = 2,
      n       = 200,
      p       = 1,
      param   = list(0, 0),
      phi_lv  = matrix(c(0.5, 0.1, 0.1, 0.4), 2, 2),
      family  = list("Gaussian", "Gaussian")
    ),
    tc_family = list("Gaussian", "Gaussian")
  )

)

# ── Run one scenario, return captured outputs ────────────────────────────────────

run_scenario <- function(sc) {

  # Stage 1: simulation
  set.seed(sc$seed)
  sim <- do.call(latent_var_sim, sc$sim_args)

  # Stage 2: timecop constructor (pipeline: observed cov → latent cov → Yule-Walker matrices)
  tc_obj <- timecop(
    data   = t(sim$X_t),
    family = sc$tc_family,
    p      = 1
  )

  # Stage 3: fit
  fit <- fit_timecop(tc_obj)

  list(
    # Simulation outputs
    sim_X_t      = sim$X_t,
    sim_Z_t      = sim$Z_t,
    sim_A_true   = sim$A_true,
    sim_Sigma    = sim$Sigma_true,

    # Constructor: covariance pipeline outputs
    cov_x_hat    = tc_obj@cov_x_hat,
    cov_z_hat    = tc_obj@cov_z_hat,
    gamma_hat    = tc_obj@gamma_hat,
    Gamma_hat    = tc_obj@Gamma_hat,

    # Fit outputs
    estimates    = fit$estimates,
    se           = fit$se
  )
}

# ── Generate mode ───────────────────────────────────────────────────────────────

generate_reference <- function() {
  cat_rule("=")
  cat("GENERATING REFERENCE SNAPSHOT\n")
  cat_rule("=")
  cat(sprintf("Output file : %s\n", REFERENCE_FILE))
  cat(sprintf("Tolerance   : %g\n\n", TOLERANCE))

  reference <- list()

  for (nm in names(SCENARIOS)) {
    sc <- SCENARIOS[[nm]]
    cat(sprintf("Running: %s ... ", sc$label))
    reference[[nm]] <- tryCatch(
      run_scenario(sc),
      error = function(e) {
        cat("ERROR\n")
        stop(sprintf("Scenario '%s' failed during generation:\n  %s", nm, e$message))
      }
    )
    cat("done\n")
  }

  saveRDS(reference, REFERENCE_FILE)
  cat(sprintf("\nReference saved to: %s\n", REFERENCE_FILE))
  cat("Commit this file to version control to share the baseline.\n")
  cat_rule()
}

# ── Compare mode ────────────────────────────────────────────────────────────────

compare_reference <- function() {
  cat_rule("=")
  cat("REGRESSION CHECK\n")
  cat_rule("=")
  cat(sprintf("Reference   : %s\n", REFERENCE_FILE))
  cat(sprintf("Tolerance   : %g\n\n", TOLERANCE))

  if (!file.exists(REFERENCE_FILE)) {
    stop(paste0(
      "Reference file not found: ", REFERENCE_FILE, "\n",
      "Run with --generate first to create it."
    ))
  }

  reference <- readRDS(REFERENCE_FILE)
  all_pass  <- TRUE

  for (nm in names(SCENARIOS)) {
    sc <- SCENARIOS[[nm]]
    cat_rule("-")
    cat(sprintf("Scenario: %s\n", sc$label))

    # Check scenario exists in reference
    if (!nm %in% names(reference)) {
      pass_fail(FALSE, sprintf("Scenario '%s' missing from reference file", nm))
      all_pass <- FALSE
      next
    }

    ref <- reference[[nm]]

    # Run current code
    current <- tryCatch(
      run_scenario(sc),
      error = function(e) {
        pass_fail(FALSE, sprintf("scenario threw an error: %s", e$message))
        NULL
      }
    )
    if (is.null(current)) { all_pass <- FALSE; next }

    # Compare each captured quantity
    checks <- list(
      list(key = "sim_X_t",    label = "latent_var_sim: X_t (observed data)"),
      list(key = "sim_Z_t",    label = "latent_var_sim: Z_t (latent data)"),
      list(key = "sim_A_true", label = "latent_var_sim: A_true (scaled coefficients)"),
      list(key = "sim_Sigma",  label = "latent_var_sim: Sigma_true (innovation covariance)"),
      list(key = "cov_x_hat",  label = "timecop:        cov_x_hat (observed covariance array)"),
      list(key = "cov_z_hat",  label = "timecop:        cov_z_hat (latent covariance array)"),
      list(key = "gamma_hat",  label = "timecop:        gamma_hat (stacked lag covariances)"),
      list(key = "Gamma_hat",  label = "timecop:        Gamma_hat (Toeplitz matrix)"),
      list(key = "estimates",  label = "fit_timecop:    estimates (VAR coefficient matrix)"),
      list(key = "se",         label = "fit_timecop:    se (standard errors)")
    )

    for (chk in checks) {
      ok <- num_equal(current[[chk$key]], ref[[chk$key]], chk$label)
      if (!ok) all_pass <- FALSE
    }
  }

  cat_rule("=")
  if (all_pass) {
    cat("\033[32mAll checks passed.\033[0m\n")
  } else {
    cat("\033[31mOne or more checks failed. See details above.\033[0m\n")
  }
  cat_rule("=")

  invisible(all_pass)
}

# ── Entry point ─────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

load_timecop()

if ("--generate" %in% args) {
  generate_reference()
} else {
  ok <- compare_reference()
  if (!ok) quit(status = 1L)
}
