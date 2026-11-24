# Multitask Code Review — Accuracy Assessment & Improvement Notes

Date: 2026-07-06
Scope: all multitask code in the package —
`timecopObjectClass.R` (multitask class + constructor), `multitask_pgd.R`,
`multitask_ebic.R`, `fit_multitask.R`, `multitask_glasso_admm.R`,
`multitask_gvar_alternate.R`, `multitask_gvar_ebic.R`,
`multitask_adaptive_weights.R`, plus the multitask tests.

Method: every estimator component was re-derived by hand and compared line by
line against the implementation, and the analytic review was backed by targeted
numerical experiments (documented inline below). Verdict up front: **no
mathematical/estimation errors were found**. The issues are one API
inconsistency, one confirmed reproducibility defect, several missing
guardrails, and a set of refactoring/performance opportunities.

---

## Part A — What was verified correct (and how)

### A1. Precision-weighted temporal objective and gradient (`multitask_pgd.R`)

The note's temporal objective for task k is (N_k/2) tr{S_eps,k(B_k) Omega_k}.
Expanding S_eps,k and dropping the B-independent term tr(C_k Omega_k), the
B-dependent part is

    N_k [ 1/2 tr(B G B' Omega) - tr(H B' Omega) ]

(the two cross terms collapse into one because Omega is symmetric:
tr(B H' Omega) = tr(H B' Omega)). The solver works in transposed space
(P = B'), with G_arg = N_k S0 and H_arg = N_k t(S10) already carrying the N_k
factor. Substituting P = B' into the code's expressions:

- smooth term: `0.5*sum((G %*% Ck) * (Ck %*% Omega)) - sum(Ck * (H %*% Omega))`
  equals N_k [ 1/2 tr(B G_nat B' Omega) - tr(H_nat B' Omega) ]  — exact match,
  verified by cyclic-trace manipulation.
- gradient: d/dP [ 1/2 tr(P' G P Omega) - tr(P' H Omega) ] = (G P - H) Omega —
  exactly the code's `(G %*% Ck - H) %*% Omega[[k]]`.

The Omega = NULL path is the identity special case and was previously verified
byte-identical to an explicit identity Omega, and to per-subject OLS at
lambda = 0 (max error ~1e-6 across tests).

### A2. Power iteration validity in the weighted case (`multitask_pgd.R`)

`apply_Q` is precisely the Hessian operator of the stacked quadratic form
f(V) = sum_k 1/2 tr((V0+Vk)' G_k (V0+Vk) Omega_k): its gradient blocks are
G_k (V0+Vk) Omega_k with the common block summed. A Hessian of a quadratic is a
symmetric PSD linear operator, so power iteration on it converges to the
correct largest eigenvalue (the Lipschitz constant) in the weighted case too.

### A3. Residual-covariance call convention (all call sites)

`gvar_resid_cov(A, S0_plus, S0_minus, S01, S10)` computes
S0 - A S01 - S10 A' + A S0 A', where S10 must be the CURRENT-PAST covariance
Cov(z_t, z_{t-1}) and S01 its transpose. In the multitask code, the local
`S10 = cov_z_hat[,,p]` is current-past, and all four call sites pass
`(B, S0, S0, t(S10), S10)`:

- `multitask_ebic.R` (pooled Sigma loop)
- `multitask_gvar_ebic.R` (joint likelihood loop)
- `multitask_gvar_alternate.R` (`resid_cov` helper)
- `fit_multitask.R` (gvar-branch grid setup)

All four agree with the required convention. Consistent.

### A4. Precision ADMM updates (`multitask_glasso_admm.R`)

- The per-task Omega update solves w(-log det O + tr(S O)) + rho/2 ||O - A||^2.
  First-order condition -w O^{-1} + w S + rho(O - A) = 0 rearranges to
  O - (w/rho) O^{-1} = A - (w/rho) S = C; eigendecomposing C and solving the
  scalar quadratic per eigenvalue gives theta_j = (c_j + sqrt(c_j^2 + 4w/rho))/2,
  all strictly positive, hence PD by construction. Matches the code exactly.
- The edgewise common-plus-unique coordinate updates
  e_k <- soft(a_k - m, lambda_E W_Ek(i,j)/rho) and
  m <- soft(mean_k(a_k - e_k), lambda_M W_M(i,j)/(rho K))
  are the exact scalar KKT solutions of the per-edge problem
  lambda_M W_M |m| + lambda_E sum_k W_Ek |e_k| + rho/2 sum_k (a_k - m - e_k)^2.
  (The /K on the m threshold arises because m enters K squared terms; the
  minimizer over m of rho/2 sum_k (r_k - m)^2 + c|m| is
  soft(mean(r_k), c/(rho K)).) Correct, including the adaptive-weight variant.
- Warm starting cannot change the fixed point (ADMM converges to the same
  optimum from any initialization for any rho > 0); verified numerically:
  warm-starting from the converged state reproduces the same solution
  (max |Delta M| ~ 5e-10) in 1 iteration vs 7 cold.

### A5. Scaling consistency across the two blocks

The temporal block uses N_k-weighted cross-products (G = N_k S0, H = N_k t(S10));
the precision block uses NORMALIZED residual covariances with likelihood
weights w_k = N_k/2. Substituting shows both blocks minimize the SAME joint
objective sum_k (N_k/2)[-log det Omega_k + tr(S_eps,k Omega_k)] + penalties, so
the alternation is coherent — no hidden scale mismatch between blocks.

### A6. Joint-objective monotonicity (numerical check, 2026-07-06)

A manual alternation (external to the package functions, recomputing the full
penalized joint objective each outer iteration on simulated K=3, d=4 data) was
non-increasing at every step and stabilized. This is the end-to-end signal that
the two conditional solvers are optimizing the same function and the coupling
(B -> S_eps -> Omega -> weighted B) is implemented consistently.

### A7. EBIC formulas

- `multitask_ebic` (gvar = FALSE): pooled Sigma = (1/N) sum_k N_k S_eps,k, and
  the profiled Gaussian log-likelihood -N/2 (log det Sigma + d); the "+ d" is
  tr(Sigma^{-1} Sigma_hat) = d at the plug-in MLE. df counts all entries of mu
  and the Delta_k; candidate counts d^2 and K d^2. Standard Foygel–Drton form.
- `multitask_gvar_ebic` (gvar = TRUE): per-task likelihood
  sum_k (N_k/2)(log det Omega_k - tr(S_eps,k Omega_k)); temporal df over all
  entries, precision df over upper-triangle off-diagonals of M and E_k.
  The E_k diagonals are free parameters not counted in df — they are constant
  across every model on the grid, so they cannot affect the argmin; harmless.

### A8. Orientation discipline

External interfaces (weights, warm starts, outputs) are all natural-orientation;
transposition is confined to the inside of `multitask_pgd`. This was pinned
earlier by asymmetric-truth tests (a true [1,2] edge recovered at [1,2], not
[2,1]) and remains intact in the current code.

---

## Part B — Issues found (ranked by importance)

### B1. API inconsistency: `obj` means two different things

**Where:** `fit_multitask.R` — gvar = FALSE results list (`obj = object`, the
data object) vs gvar = TRUE results list (`obj = gfit$obj`, the scalar
penalized objective value, with the data object under `object`).

**Why it matters:** the same field name in the two return classes holds
completely different types. Downstream code written against one branch
(e.g. `fit$obj@subjects`) breaks silently or confusingly on the other. It also
diverges from the package's own convention: `timecop_fit` and `timecop_gvar`
both use `obj` for the embedded data object.

**Detail:** this arose historically — the gvar branch added `obj` for the
objective value (mirroring `multitask_gvar_alternate`'s internal return) and
`object` for the data, while the earlier gvar = FALSE branch already used `obj`
for the data.

**Proposed fix (small, breaking-but-pre-release):** in the gvar = TRUE result,
rename the scalar to `objective` and the data object to `obj`, matching the
other three fit classes. Update the roxygen `@return` accordingly and the one
test that touches these names (none currently reference `obj`/`objective` in
the gvar fit, so the change is low-risk).

### B2. Confirmed reproducibility defect: fits consume the user's RNG stream

**Where:** `multitask_pgd.R` — the Lipschitz power iteration starts from
`stats::rnorm((K+1) d^2)`.

**Why it matters:** every call to `fit_multitask` (both branches) advances the
global RNG state. Confirmed numerically (2026-07-06): with `set.seed(123)`,
`rnorm(1)` returns a different value depending on whether a fit ran in between.
In simulation studies — exactly this project's workflow — inserting or removing
a fit call silently changes all subsequently generated data, breaking exact
reproducibility of scripts and making "same seed, same results" false.

**Detail:** the random start is only used to seed power iteration; the fitted
estimates do not depend on it (the problem is convex and L only sets the step
size), but the RNG side effect is real, and in principle a slightly different
L can change iteration counts/paths.

**Proposed fix:** replace the random start with a deterministic full-spectrum
vector (e.g. `V <- matrix(sin(seq_len((K+1)*d*d)), ...)`, normalized). A fixed
generic vector almost surely has a component along the top eigenvector, which
is all power iteration needs, and it makes the solver a pure function.
(Alternative: save/restore `.Random.seed`, but the deterministic vector is
simpler and also makes L itself reproducible.)

### B3. No safety margin on the Lipschitz estimate + restart accepts uphill steps

**Where:** `multitask_pgd.R` — power iteration result used as-is (`1/L` step);
in the monotone-restart branch, the ISTA fallback's `Bt_new` is accepted
unconditionally, even if its objective is still higher than `obj_old`.

**Why it matters:** power iteration converges to lambda_max FROM BELOW; if it
stops early (tolerance 1e-6, cap 100 iterations), L slightly underestimates the
true constant and the gradient step 1/L is slightly too large. The FISTA
restart guards the accelerated step, but the fallback ISTA step is then taken
on faith. With an underestimated L, both steps can overshoot, and the loop
would accept an objective increase — in pathological cases this can oscillate
rather than converge. Never observed in our tests, but it is a hole.

**Proposed fix:** multiply the estimate by a small safety factor
(`L <- 1.01 * L`) after power iteration. One line, negligible cost, closes the
overshoot risk. (A stricter alternative — reject the fallback step if it
increases the objective and halve the step — is more code for a case the
safety factor already covers.)

### B4. `penalty = "scad"` with `gvar = TRUE` silently degrades to lasso

**Where:** `fit_multitask.R` — the weight block only computes weights for
`penalty == "adaptive"`; the SCAD LLA loop exists only in the gvar = FALSE
`fit_point()`. In the gvar branch, `penalty = "scad"` leaves all `W_*` NULL, so
the alternation runs a plain (uniform) lasso while the result object reports
`penalty = "scad"`.

**Why it matters:** the user asks for SCAD and gets lasso with no signal —
worse, the returned object *claims* SCAD. This is a silent correctness lie in
the output metadata.

**Proposed fix:** until/unless SCAD is wired into the alternation (LLA around
the outer loop), `stop()` with a clear message ("SCAD is not yet supported with
gvar = TRUE; use 'lasso' or 'adaptive'"). An error is better than a warning
here because the output metadata would otherwise misreport the penalty.

### B5. Non-convergence is silent everywhere in the gvar path

**Where:** `multitask_glasso_admm` returns a `converged` flag that
`multitask_gvar_alternate` never inspects; the outer alternation hitting
`max_outer` without meeting `tol` is also silent; `multitask_pgd` hitting
`max_iter` likewise.

**Why it matters:** a grid point that never converged contributes its EBIC to
the selection as if it were a converged fit. If the *selected* model happens to
be one of these, the user gets no indication the estimates are unpolished.

**Proposed fix:** track convergence through the alternation (e.g. return
`converged = admm_converged && outer_converged`) and issue ONE `warning()` in
`fit_multitask` if the EBIC-selected best fit did not converge. Avoid warning
per grid point (noise).

### B6. `lambda_M` / `lambda_E` cannot be vectors (inconsistent with temporal penalties)

**Where:** `fit_multitask.R` validation (`length(lambda_M) != 1` rejected)
vs. `lambda_mu`/`lambda_delta`, which accept vectors as custom grids.

**Why it matters:** the 4-D grid code itself already handles vector
`lambda_M`/`lambda_E` (`sort(lambda_M, decreasing = TRUE)` becomes the axis);
only the validator forbids it. Users who want a custom precision grid (e.g. a
finer sweep near a known region) can do it for the temporal axes but not the
precision axes, for no reason.

**Proposed fix:** relax the two validators to "non-negative numeric scalar or
vector" (same wording as the temporal ones), and update the `@param` docs.

### B7. Wasted / unguarded computation and inconsistent ridging in the gvar branch

**Where:** `fit_multitask.R` gvar branch: `B_yw0`, `S_eps0`, `Om0` are computed
unconditionally, even when both `lambda_M` and `lambda_E` are supplied (their
only purpose is the data-driven precision grid). Also `B_yw0` uses plain
`solve(S0[[k]])` while `Om0` uses `solve(S_eps0 + 1e-3 I)`; the unpenalized
pilot in `multitask_adaptive_weights` likewise uses plain `solve(S0[[k]])`.

**Why it matters:** (i) minor wasted work; (ii) more importantly, if a
subject's latent lag-0 matrix is near-singular (short series, strong
discreteness, no `pd_approx`), the plain `solve()` errors out even when the
user supplied both precision penalties and the quantity was never needed;
(iii) the inconsistent use of ridge regularization (some solves guarded, some
not) means failure modes differ by code path for the same data.

**Proposed fix:** wrap the grid-derivation block in
`if (is.null(lambda_M) || is.null(lambda_E))`, and use a single ridged solve
helper (e.g. `solve(S + 1e-8 * mean(diag(S)) * I)`) consistently for all
pilot/YW solves in the multitask code.

### B8. Stale documentation from the removed staged selector; grid-max nuance

**Where:** `fit_multitask.R` roxygen — `@param lambda_M` / `@param lambda_E`
say "selected by EBIC (Stage B)"; `@param n_lambda_prec` says "Stage-B ...
EBIC search". The staged selector was replaced by the full 4-D grid; there are
no stages anymore.

Additionally (nuance, not an error): the temporal grid tops
`lam_mu_max = max|sum_k H_k|` and `lam_d_max = max_k max|H_k|` are the exact
zeroing thresholds for the UNWEIGHTED (Omega = I) problem. Under precision
weighting the gradient at zero is `-H Omega`, so the exact full-sparsity
threshold shifts; the largest grid value may not fully zero the temporal
networks in the gvar case. Selection still works (the grid brackets the
interesting region), but the top-of-grid interpretation differs slightly
between branches.

**Proposed fix:** reword the three `@param`s to describe the 4-D grid; add one
sentence to the gvar `@param` noting the grid anchors are derived from the
unweighted problem.

### B9. No finiteness guard in the EBIC functions

**Where:** `multitask_ebic.R` and `multitask_gvar_ebic.R` — `determinant()` of
the pooled/implied covariance and the estimated precisions.

**Why it matters:** `S_eps,k(B)` is PSD only if the joint latent covariance
block matrix is PSD; interpolated latent covariances can violate this (the
package clamps to [-1,1] and offers `pd_approx`, but neither guarantees the
lag-0/lag-1 block structure is jointly PSD). A non-PD pooled Sigma makes
`determinant()$modulus` NaN; the NaN EBIC then propagates through the grid
comparison (`NaN < best` is NA -> `if()` error, or silently never selected
depending on position). The failure mode is confusing for the user.

**Proposed fix:** after computing each EBIC, `if (!is.finite(ebic)) ebic <- Inf`
(with an optional single message). Inf cleanly loses every grid comparison and
turns a cryptic failure into "that grid point is excluded".

---

## Part C — Refactoring & performance improvements

### C1. Vectorize the edgewise ADMM decomposition (largest remaining perf lever)

**Where:** `multitask_common_unique_update` (`multitask_glasso_admm.R`) — a
double `for` over d(d-1)/2 edges, each with an inner coordinate-descent loop,
in pure R, executed EVERY ADMM iteration, which itself runs every outer
alternation of every grid point.

**Why it matters:** the per-edge problems are independent and their updates are
elementwise operations, so the entire sweep can be done on whole matrices:

    repeat {
      E_k <- soft(R_k - M, T_E_k)          # elementwise, all edges at once
      M   <- soft_offdiag(mean_k(R_k - E_k), T_M / K)
    } until max change < tol

with T_E_k / T_M the (possibly weighted) threshold matrices. This is the SAME
fixed point (it is exactly the same alternation, applied to all edges
simultaneously instead of one at a time), but removes the O(d^2 K) interpreter
overhead per pass. For d = 4 the current loop is fine; for d = 15–30 (realistic
EMA panels) the pure-R edge loop will dominate the entire fit.

**Care needed:** keep diag(M) = 0 and the E_k-diagonal absorption exactly as
now; verify equivalence against the current implementation on random problems
(same M/E to ~1e-8) before swapping in.

### C2. Extract the gvar branch out of `fit_multitask`

**Where:** `fit_multitask.R` is ~445 lines; the `if (gvar) { ... }` block is
~100 lines with its own grid, warm-start state, and result assembly.

**Why it matters:** readability and testability. The method currently does
validation + block extraction + two entirely different selection algorithms.
Moving the 4-D grid into an internal `multitask_gvar_grid(...)` (same file or
its own) makes each piece unit-testable and the S4 method a thin dispatcher.
Pure code motion; zero behavior change.

### C3. DRY and micro-efficiency inside `multitask_pgd`

- `apply_Q(V)` is exactly `grad_stacked(V)` with H = 0. Unify into one
  function with an optional `include_H` flag (or have grad call the operator
  and subtract the constant term). Removes ~15 duplicated lines and one
  maintenance hazard (the two must stay in sync w.r.t. the Omega weighting —
  they already diverged once during development).
- `tau / L` is recomputed (a full matrix allocation) at every FISTA iteration;
  hoist `tau_L <- tau / L` before the loop.
- When `Omega = NULL` and G is fixed (the entire gvar = FALSE grid), L is
  identical across all grid points but recomputed per call. Optionally accept a
  precomputed `L` argument and compute it once in `fit_multitask`. Minor.

### C4. Expose `n_lambda_init` for the penalized pilot

**Where:** `multitask_adaptive_weights(..., n_lambda_init = 10L)` is not
reachable from `fit_multitask`.

**Why it matters:** the penalized pilot runs a 10x10 EBIC grid per subject
(K fits). The pilot only needs a rough consistent estimate; users with many
subjects may want `n_lambda_init = 4` (or a single lambda). One passthrough
argument.

### C5. Add `print()` / `summary()` methods for the two multitask fit classes

**Where:** `timecop_multitask_fit` and `timecop_multitask_gvar_fit` have no
methods; printing dumps the entire list INCLUDING the embedded
`timecop_multitask` object (all K subjects' data matrices) — pages of noise.
The single-subject classes (`timecop_fit`, `timecop_gvar`) have proper methods
in `R/summary.R`.

**Proposed:** compact `print` methods (dimensions, K, selected penalties, df,
EBIC, convergence) and `summary` methods mirroring the existing style —
e.g. common network edge table, per-subject deviation counts.

### C6. Test-coverage gaps

Current `test-multitask.R` (28 assertions) covers the constructor, PGD OLS
recovery, Omega = NULL equivalence, ADMM basics, both fit branches, and grid
collapse. Not covered:

- adaptive gvar path (`penalty = "adaptive"` with gvar = TRUE, both
  `adaptive_weights_from` values) — currently only exercised manually;
- ADMM warm-start-equals-cold-solution equivalence (regression-guard for the
  warm-start code);
- `multitask_gvar_ebic` df counting on a hand-constructed case (e.g. known
  nonzero patterns -> exact df vector);
- vectorized-vs-edgewise equivalence if C1 is implemented;
- an RNG-purity test if B2 is fixed (`.Random.seed` unchanged by a fit).

---

## Suggested implementation order

1. **B2, B3** — one-file fixes in `multitask_pgd.R` (deterministic power-iteration
   start + L safety factor). Lowest risk, immediate reproducibility win.
2. **B1** — rename `obj`/`objective` in the gvar result (API consistency).
3. **B4, B5, B9** — guardrails (SCAD error, non-convergence warning, finite-EBIC).
4. **B6, B7, B8** — validator relaxation, guarded/ridged solves, doc cleanup.
5. **C1** — vectorized decomposition (with an equivalence test), then C3.
6. **C2, C4, C5, C6** — structural extraction, pilot knob, print methods, tests.
