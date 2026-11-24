# Diagnosis: Dense Recovered Precision Matrices in the Multitask Graphical VAR

Date: 2026-07-06
Companion to: `multitask_code_review.md`
Status: diagnosis confirmed numerically; **fix APPLIED and validated 2026-07-06**.
Implementation: grid construction centralized in `R/multitask_lambda_grids.R`
(`multitask_prec_anchors()` = the Sec. 4 KKT math in one place;
`multitask_lambda_grids()` = all four sequences), wired into `fit_multitask`
and `multitask_gvar_alternate`'s fallback defaults; boundary-of-auto-grid
selections now emit a message (both branches); KKT anchor test added to
`test-multitask.R` (grid tops must produce fully sparse fits).
Validation result (one-true-edge scenario, Sec. 3 setup, full API): lasso
selects interior lambda_M = 90.5 and returns nnz_M = 1 (the true edge,
M[1,3] = -0.36 vs true -0.35), nnz_E = 0; adaptive returns df = (5, 0, 1, 0) —
exact temporal and precision supports; temporal recovery unchanged; full test
suite passes (33 multitask assertions).

---

## 1. Symptom

`fit_multitask(gvar = TRUE)` recovers the temporal networks (mu, Delta_k) well,
but the recovered precision networks (M, E_k, hence Omega_k) are consistently
**too dense** — including under `penalty = "adaptive"` with unpenalized pilot
weights. A telltale signature was visible in every earlier test run: the
EBIC-selected precision penalties always landed at or near the very BOTTOM of
their grid (e.g. lambda_M = 0.00062, 0.00066, lambda_E = 0.0025–0.0066), while
the selected temporal penalties landed at interior grid points.

Three candidate culprits were considered:

1. the lambda grid (are the candidate penalties in the right range?),
2. the optimization (does the ADMM solve its subproblem correctly?),
3. the EBIC (does the criterion over-favor dense precision models?).

## 2. Verdict (one paragraph)

**The lambda grid is the root cause.** The precision-penalty grid is anchored
about two orders of magnitude below the range in which the penalty has any
sparsifying effect, so every candidate model on the grid is dense; EBIC can
only choose among the candidates it is offered, and it dutifully selects the
least-bad dense one (the highest-likelihood fit, at the bottom of the grid).
The ADMM optimizer is correct — evaluated at a properly scaled penalty it
recovers exactly the true precision support. The EBIC is (mostly) innocent,
with two second-order caveats noted in Section 8.

## 3. The diagnostic experiment

Setup: K = 3 Gaussian subjects, d = 4, n = 800 each, common temporal matrix
(diag 0.4, one off-diagonal edge), and exactly ONE true shared precision edge
(Omega = I + M_true with M_true[1,3] = M_true[3,1] = -0.35). Per-subject
residual covariances S_eps,k were formed from the Yule-Walker transitions, and
`multitask_glasso_admm` was run at penalties spanning both the coded grid range
and the theoretically derived range, with likelihood weights w_k = N_k / 2
(as in the actual pipeline).

Measured quantities:

    grid anchor as coded  : 0.5647     (= max |offdiag(Om0)|, see Sec. 5)
    grid bottom as coded  : 0.000565   (= anchor * 1e-3)
    KKT anchor for E      : 126.8      (derived scale, see Sec. 4)
    KKT anchor for M      : 360.4

Sparsity of the ADMM solution as lambda (= lambda_M = lambda_E) increases,
counting nonzero upper-triangle off-diagonals (true values: nnz_M = 1,
nnz_E = 0):

    lambda =   0.001   ->  nnz_M = 4,  nnz_E = 17    (dense)
    lambda =   0.565   ->  nnz_M = 4,  nnz_E = 18    (dense)  <- TOP of coded grid
    lambda =   6.3     ->  nnz_M = 4,  nnz_E = 11
    lambda =  25       ->  nnz_M = 2,  nnz_E = 0
    lambda = 127       ->  nnz_M = 1,  nnz_E = 0     <- exactly the true support

Interpretation: across the ENTIRE coded grid [0.00057, 0.565] the solution
barely changes and is dense everywhere — the penalty is effectively zero on
that range. All sparsification happens between lambda ~ 6 and ~ 127, i.e. one
to two orders of magnitude ABOVE the top of the coded grid. At the derived KKT
anchor the solver returns exactly the true single shared edge, which
simultaneously certifies that the ADMM optimization is correct.

## 4. Why: the KKT zeroing scale carries the likelihood weight N_k/2

The precision block minimizes

    sum_k w_k [ -log det(Omega_k) + tr(S_eps,k Omega_k) ]
      + lambda_M ||M||_{1,off} + lambda_E sum_k ||E_k||_{1,off},
    Omega_k = M + E_k,   w_k = N_k / 2.

Consider when the fully sparse solution (all off-diagonals of M and E_k zero,
i.e. diagonal Omega_k) satisfies the first-order conditions. The gradient of
the smooth part with respect to the off-diagonal entry (Omega_k)_ij is

    w_k ( S_eps,k - Omega_k^{-1} )_ij     (per symmetric-pair entry; the
                                           factor-2 from double counting the
                                           (i,j)/(j,i) pair appears on both the
                                           loss and the ||.||_{1,off} penalty
                                           and cancels).

At a diagonal Omega_k, (Omega_k^{-1})_ij = 0 off the diagonal, so the
subgradient condition for keeping edge (i,j) of E_k at zero is

    | w_k * (S_eps,k)_ij |  <=  lambda_E                       (E-edge)

and for the shared entry m = M_ij, whose gradient sums across tasks,

    | sum_k w_k * (S_eps,k)_ij |  <=  lambda_M.                (M-edge)

Therefore the smallest penalties that produce a fully sparse precision network
— the natural TOPS of the two grids, exactly analogous to
`lambda_max = max|X'y|` in lasso regression — are

    lam_E_max = max_k [ w_k * max_offdiag |S_eps,k| ]
    lam_M_max = max_offdiag | sum_k w_k * S_eps,k |

Both scales are proportional to N_k/2 (hundreds, for typical series lengths),
because the likelihood term the penalty must overcome is itself weighted by
N_k/2. In the experiment: lam_E_max ~ 127, lam_M_max ~ 360 — matching where
sparsification was empirically observed.

## 5. What the code anchors on instead (the two compounding errors)

`fit_multitask.R` (gvar branch) builds ONE shared anchor for both precision
grids:

    Om0          <- solve(S_eps0 + 1e-3 I)             # initial precision
    lam_prec_max <- max |offdiag(Om0)|                  # ~ 0.1 - 1
    prec_seq     <- logspace(lam_prec_max, lam_prec_max * 1e-3, n_lambda_prec)

Two independent errors compound:

1. **Wrong matrix.** The anchor uses the off-diagonals of the estimated
   PRECISION matrix (entries ~ 0.1-0.6, i.e. partial-covariance scale), whereas
   the KKT threshold involves the off-diagonals of the residual COVARIANCE
   S_eps,k. (In single-task graphical lasso the standard anchor is
   max|offdiag(S)| — of the covariance — precisely because of the same KKT
   argument.)
2. **Missing likelihood weight.** Even with the right matrix, the coded anchor
   omits the w_k = N_k/2 factor. The single-task convention it was ported from
   (`fit_graphical_var` -> glasso) uses an UNWEIGHTED likelihood
   (-log det + tr(S Theta)), where an S-scale anchor is correct. Our multitask
   objective multiplies the likelihood by N_k/2, so the penalty must scale up
   by the same factor to compete. The N never came along in the port.

Combined effect in the experiment: anchor 0.56 vs. required ~127-360, a factor
of ~200-600. The log-spaced grid descends a further 1e-3 from its anchor, so
ALL candidate penalties sit deep inside the no-effect region.

## 6. Why the temporal side does NOT have this problem

This is the key to the observed asymmetry ("temporal effects recover well,
precision matrices don't"). The temporal grids are anchored on

    lam_mu_max    = max | sum_k H_k |,      H_k = N_k * t(S10_k)
    lam_delta_max = max_k max | H_k |

and H_k already CARRIES the N_k factor (it is the N-weighted cross-product fed
to the PGD solver). These are exactly the KKT zeroing thresholds for mu and
Delta_k in the unweighted-precision (Omega = I) problem: the gradient of the
temporal loss at zero coefficients is -H_k (blockwise), so full temporal
sparsity occurs at max|sum_k H_k| for mu — precisely the coded anchor. In other
words, the temporal grid was built on the correctly weighted quantity from the
start, while the precision grid was not. Same estimator, one calibrated axis
pair, one miscalibrated axis pair — hence good temporal recovery and dense
precision recovery from the same fit.

## 7. Why adaptive weights sometimes masked the problem

Adaptive weights multiply each edge's effective penalty by
W = 1/(|pilot| + eps)^gamma. For a null edge whose pilot estimate is
noise-level (say |pilot| ~ 0.01-0.001), W ~ 100-1000, which can accidentally
lift lambda * W into the sparsifying range EVEN THOUGH lambda itself is far too
small. This is why some earlier small-d sanity tests with penalty = "adaptive"
recovered df_M = 1 and looked healthy: the weights performed an uncontrolled,
data-dependent rescue of a miscalibrated grid. The rescue is unreliable —
it depends on how small the pilot's null entries happen to be — which is
consistent with the user's finding that adaptive + unpenalized pilot still
yields overly dense precision matrices on other data. Calibration must come
from the grid; adaptive weights should only provide RELATIVE (per-edge)
discrimination around it.

## 8. The EBIC's role (mostly innocent, two second-order caveats)

EBIC selects the minimizer over the offered grid. When every offered model is
dense, EBIC's choice is confined to dense models, and among them the likelihood
term dominates (adding N_k/2-weighted likelihood for near-MLE fits), pushing
the selection to the smallest offered penalty — exactly the observed
bottom-of-grid selections. Given a correctly ranged grid (Section 9), the same
EBIC selected sparse, correct supports in our verification runs, so the
criterion itself is not the driver.

Two genuine but second-order EBIC caveats to revisit AFTER the grid fix:

- **df counting at tol_zero = 1e-8.** Any noise-level nonzero counts as an
  edge. With a properly scaled penalty most such entries are exactly zero (soft
  thresholding), so this mostly self-resolves; but at small penalties the df is
  inflated and, more subtly, an edge can be double-counted when its value
  splits between m and one or more e_k (df = 1 + #tasks for what is
  functionally one edge). If over-selection persists post-fix, consider
  counting an edge once per task via nnz(Omega_k) or raising tol_zero.
- **Penalized-likelihood plug-in.** The EBIC evaluates the likelihood at the
  shrunken (penalized) estimates rather than at refit MLEs on the selected
  support. This biases comparisons slightly toward smaller penalties for both
  blocks; standard practice, but a refit-then-score variant would remove it if
  needed.

## 9. Proposed fix

Anchor the two precision grids at their KKT scales, with SEPARATE sequences for
M and E (their scales genuinely differ: M's anchor aggregates over tasks):

    w         <- N / 2
    S_eps0_k  <- gvar_resid_cov(B_yw_k, S0_k, S0_k, t(S10_k), S10_k)
    lam_E_max <- max_k ( w_k * max|offdiag(S_eps0_k)| )
    lam_M_max <- max|offdiag( sum_k w_k * S_eps0_k )|
    lambda_E_seq <- logspace(lam_E_max, lam_E_max * 1e-3, n_lambda_prec)
    lambda_M_seq <- logspace(lam_M_max, lam_M_max * 1e-3, n_lambda_prec)

Notes on the fix:

- This mirrors exactly how the temporal anchors already work, restoring
  symmetry between the two blocks.
- The same correction applies to the fallback defaults inside
  `multitask_gvar_alternate` (currently 0.1 * max|offdiag(S_eps)| WITHOUT the
  w_k factor — same missing weight).
- The anchors should be computed from the Yule-Walker initial S_eps (as now),
  which is the standard practice (anchor from a pilot; the exact anchor shifts
  slightly as B_k updates during the alternation, which is immaterial for a
  log-spaced grid).
- Grid depth: descending 1e-3 from a now-correct anchor spans from full
  sparsity down to effectively unpenalized — the full path. If grids feel
  wasteful at the bottom, 1e-2 is a reasonable floor.
- Under adaptive weights the effective per-edge threshold is lambda * W; no
  change to the anchor is needed (weights are centered near 1 for
  moderate pilot entries and provide relative discrimination around the
  calibrated lambda).

## 9b. Follow-up (2026-07-07): weight-aware anchors for adaptive penalties

After the Section-9 fix, ADAPTIVE fits still selected lambda_M at the BOTTOM of
the auto grid (dense end, df_M inflated by one false positive). Cause: the
anchors implement the UNIFORM-lasso KKT condition |g_ij| <= lambda, but under
adaptive weights the condition is |g_ij| <= lambda * W_ij, i.e. each edge zeros
at t_ij = g_ij / W_ij. Since all precision weights exceed 1 (pilot precision
entries < 1, W = 1/(|pilot|+eps)^gamma; true edges W ~ 3, nulls W ~ 50-1000),
the weighted problem's active range sits 1-3 decades BELOW the unweighted
anchor and is STRETCHED (spread of thresholds ~ spread(g) x spread(W)), so the
3-decade auto descent truncated the null-entry zone at the grid bottom. Under
lasso (W = 1) the anchors were exact, which is why lasso selected interior
points while adaptive rode the boundary.

Fix: `multitask_prec_anchors` and `multitask_lambda_grids` now accept the
weight matrices and divide the gradients ELEMENTWISE by the weights before
taking maxima (all four axes; temporal division uses t(W) because H lives in
the solver's transposed space while weights are natural-orientation).
`fit_multitask` computes the adaptive weights BEFORE building the grids and
passes them in; `multitask_gvar_alternate`'s fallback defaults are likewise
weight-aware. W = NULL reduces to the uniform anchors (lasso unchanged).

Validated: the previously boundary-riding adaptive scenario now selects
lambda_M at grid index 6/10 and lambda_E at 4/10 (both interior), with exact
structural recovery df = (5, 0, 1, 2): shared edge [1,3] alone in M, each
unique edge in the correct subject's E_k, no-deviation subject empty. A
weighted KKT anchor test was added (weighted grid tops must fully sparsify the
weighted fits; weighted anchors must sit below uniform ones).

## 10. Validation plan for the fix

1. Re-run the Section 3 experiment through the full
   `fit_multitask(gvar = TRUE)` API: EBIC should now select an interior
   (lambda_M, lambda_E) and return nnz_M = 1, nnz_E = 0 (or near), instead of
   bottom-of-grid dense fits.
2. Confirm the temporal recovery is unchanged (it should be: temporal anchors
   untouched).
3. Confirm the selected penalties land at interior grid points on repeated
   simulations (a selection at either grid END is the standard sign of a
   misranged grid and should be treated as a red flag going forward — worth a
   runtime message when it happens).
4. Re-run the adaptive variants (unpenalized and penalized pilots): both should
   now produce sparse precision networks WITHOUT relying on the accidental
   weight-rescue effect.
5. Run the existing test suite (grid-collapse and best==min(grid) tests are
   unaffected in logic but re-verify) and update any test whose supplied
   lambda_M/lambda_E values were chosen on the old scale.
