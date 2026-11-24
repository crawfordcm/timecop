# The Multitask Graphical VAR in `timecop`: Model, Estimation, and Optimization

This document is a self-contained mathematical description of the multitask
graphical VAR functionality in the `timecop` package, i.e. what happens when you
call

```r
obj <- timecop_multitask(data = list_of_matrices, family = family)
fit <- fit_multitask(obj, gvar = TRUE)
```

It covers the latent-Gaussian (copula) model for discrete multivariate time
series, the common-plus-unique decomposition of both the temporal and the
contemporaneous (precision) networks across subjects, the penalized
pseudo-likelihood that is minimized, the two solvers (an accelerated proximal
gradient method for the temporal block and an ADMM for the precision block), the
outer alternation that couples them, the data-driven construction of the penalty
grids, adaptive weights, and the joint EBIC used for model selection. Section
numbers in the margin of the code comments ("note section 5", etc.) refer to
the design notes this implementation follows; here everything is derived from
scratch.

The intended reader is a graduate student in statistics who is comfortable with
multivariate Gaussian likelihoods, convex optimization, and the lasso, but has
not seen this particular model before.

---

## Contents

1. [Notation and setting](#1-notation-and-setting)
2. [The single-subject latent Gaussian VAR](#2-the-single-subject-latent-gaussian-var)
3. [The multitask graphical VAR model](#3-the-multitask-graphical-var-model)
4. [From data to the covariance-only objective](#4-from-data-to-the-covariance-only-objective)
5. [The penalized objective](#5-the-penalized-objective)
6. [Optimization I: the outer alternation](#6-optimization-i-the-outer-alternation)
7. [Optimization II: the temporal block (FISTA)](#7-optimization-ii-the-temporal-block-fista)
8. [Optimization III: the precision block (ADMM)](#8-optimization-iii-the-precision-block-admm)
9. [Penalty grids and KKT anchors](#9-penalty-grids-and-kkt-anchors)
10. [Adaptive weights](#10-adaptive-weights)
11. [Model selection: the joint EBIC](#11-model-selection-the-joint-ebic)
12. [Searching the four-dimensional penalty lattice](#12-searching-the-four-dimensional-penalty-lattice)
13. [Special cases and relatives](#13-special-cases-and-relatives)
14. [The complete algorithm](#14-the-complete-algorithm)
15. [Map from mathematics to code](#15-map-from-mathematics-to-code)
16. [References](#16-references)

---

## 1. Notation and setting

We observe $K$ subjects (in multi-task-learning language, $K$ *tasks*). Subject
$k \in \{1,\dots,K\}$ contributes a $d$-variate time series

$$
X^{(k)}_t = \big(X^{(k)}_{t,1}, \dots, X^{(k)}_{t,d}\big)^\top, \qquad t = 1, \dots, n_k .
$$

Every subject has the same $d$ variables, and each variable $i$ has the same
marginal family across subjects (Bernoulli, Poisson, or Gaussian), but series
lengths $n_k$ may differ. Because the model is a VAR of order $p = 1$, each
subject supplies $N_k = n_k - 1$ usable (current, lagged) pairs. We write
$N = \sum_k N_k$.

Matrix conventions used throughout:

| Symbol | Meaning |
|---|---|
| $\|A\|_1 = \sum_{ij}\lvert A_{ij}\rvert$ | entrywise $\ell_1$ norm (all entries) |
| $\|A\|_{1,\mathrm{off}} = \sum_{i \ne j}\lvert A_{ij}\rvert$ | off-diagonal $\ell_1$ norm |
| $\|A\|_F$ | Frobenius norm |
| $\langle A, B\rangle = \operatorname{tr}(A^\top B) = \sum_{ij} A_{ij}B_{ij}$ | Frobenius inner product |
| $A \succ 0$ | symmetric positive definite |
| $\mathcal{S}_\tau(x) = \operatorname{sign}(x)\max(\lvert x\rvert - \tau, 0)$ | soft-thresholding, applied entrywise; $\tau$ may be a matrix of the same shape |

---

## 2. The single-subject latent Gaussian VAR

The multitask machinery is built on top of the single-subject model that the
rest of `timecop` estimates, so we describe that first. Fix one subject and drop
the superscript $(k)$.

### 2.1 Latent process

Assume there is an unobserved, stationary, zero-mean Gaussian VAR(1) process

$$
Z_t = B\, Z_{t-1} + \varepsilon_t, \qquad \varepsilon_t \overset{\text{iid}}{\sim} \mathcal{N}_d(0, \Sigma), \qquad \Omega = \Sigma^{-1},
$$

where $B \in \mathbb{R}^{d\times d}$ is the *temporal* (transition) matrix and
$\Omega$ is the *innovation precision* matrix. Stationarity requires the
spectral radius of $B$ to be below one. The latent process is normalized so
that each coordinate has unit variance, $\operatorname{Var}(Z_{t,i}) = 1$; this
is a scale convention, not a restriction (the discrete marginals below are
invariant to the scale of $Z$).

Two networks live in this model:

* the **temporal network**, whose edges are the nonzero entries of $B$: $B_{ij} \neq 0$ means variable $j$ at time $t-1$ predicts variable $i$ at time $t$;
* the **contemporaneous network**, whose edges are the nonzero off-diagonal entries of $\Omega$: $\Omega_{ij} \neq 0$ means innovations $i$ and $j$ are conditionally dependent given all other innovations at the same time (a Gaussian graphical model on the innovations).

This pairing of a sparse transition matrix with a sparse innovation precision
is the *graphical VAR*.

### 2.2 Observation model (Gaussian copula)

Each observed coordinate is a deterministic, monotone transform of the
corresponding latent coordinate,

$$
X_{t,i} = G_i(Z_{t,i}), \qquad G_i(z) = F_i^{-1}\!\big(\Phi(z)\big),
$$

where $\Phi$ is the standard normal CDF and $F_i^{-1}$ is the quantile function
of the marginal family of variable $i$ (with the convention
$F_i^{-1}(u) = \inf\{x : F_i(x) \ge u\}$ for discrete families). Concretely:

* **Bernoulli$(p_i)$**: $X_{t,i} = \mathbf{1}\{Z_{t,i} > q_i\}$ with $q_i = \Phi^{-1}(1 - p_i)$.
* **Poisson$(\lambda_i)$**: $X_{t,i} = j$ iff $\Phi^{-1}(F_i(j-1)) < Z_{t,i} \le \Phi^{-1}(F_i(j))$, i.e. the real line is cut at the normal quantiles of the Poisson CDF.
* **Gaussian**: $X_{t,i} = Z_{t,i}$. The package treats this link as the identity, so Gaussian variables **must be standardized by the user** to unit variance before they are passed in.

Marginal parameters ($p_i$ or $\lambda_i$) are estimated by the sample
proportion or sample mean of the observed series. This is a Gaussian copula
model with a VAR(1) dependence structure; see Jia, Kang and Pipiras (2023) for
the count-series version that `timecop` follows.

### 2.3 The link function: observed covariances from latent correlations

Because $X_{t,i}$ is a function of $Z_{t,i}$ alone, the observed
cross-covariance at lag $h$ between two variables depends on the latent
processes only through the latent correlation
$u = \operatorname{Corr}(Z_{t,i}, Z_{t-h,j})$. Expanding each transform in
probabilists' Hermite polynomials $\mathrm{He}_m$ (with $\mathrm{He}_0 = 1$,
$\mathrm{He}_1(z) = z$, $\mathrm{He}_2(z) = z^2 - 1$, ...),

$$
G_i(z) = \sum_{m \ge 0} g_{i,m}\, \mathrm{He}_m(z), \qquad g_{i,m} = \frac{1}{m!}\, \mathbb{E}\big[G_i(Z)\,\mathrm{He}_m(Z)\big], \quad Z \sim \mathcal{N}(0,1),
$$

and using the orthogonality relation
$\mathbb{E}[\mathrm{He}_m(Z_1)\mathrm{He}_{m'}(Z_2)] = m!\, u^m \mathbf{1}\{m = m'\}$
for a standard bivariate normal pair with correlation $u$, one obtains the
**link function**

$$
L_{ij}(u) := \operatorname{Cov}(X_{t,i}, X_{t-h,j}) = \sum_{m \ge 1} m!\, g_{i,m}\, g_{j,m}\, u^m .
$$

The coefficients $\ell_{ij,m} = m!\, g_{i,m} g_{j,m}$ have closed forms:

* **Bernoulli**: $g_{i,m} = \dfrac{\varphi(q_i)\,\mathrm{He}_{m-1}(q_i)}{m!}$ for $m \ge 1$, where $\varphi$ is the standard normal density.
* **Poisson**: $g_{i,m} = \dfrac{1}{m!}\displaystyle\sum_{j \ge 0} \varphi(q_{ij})\,\mathrm{He}_{m-1}(q_{ij})$ with $q_{ij} = \Phi^{-1}(F_i(j))$; terms with $F_i(j) \in \{0, 1\}$ are dropped, and the sum is truncated at $j = 50$.
* **Gaussian**: $g_{i,1} = 1$ and $g_{i,m} = 0$ for $m \ge 2$, so $L_{ij}(u) = u$ when both variables are Gaussian.

The package truncates the series at $m = 100$ terms and uses pre-computed
Hermite polynomials. $L_{ij}$ is a strictly increasing function on $[-1, 1]$:
by Hoeffding's covariance identity, the covariance of two nondecreasing
functions of a bivariate normal pair is nondecreasing in the correlation, and
the $m = 1$ term makes the increase strict. It is therefore invertible. When the constructor is called with
`corr = TRUE`, $L_{ij}$ is divided by the product of the marginal standard
deviations so that it maps latent correlations to observed *correlations*.

### 2.4 Inverse link and the latent covariance estimates

For $h \in \{-1, 0, 1\}$, the sample cross-covariances (or correlations) of the
observed series are computed,

$$
\hat\gamma^{X}_{ij}(h) = \frac{1}{n - \lvert h\rvert}\sum_{t} \big(X_{t,i} - \bar X_i\big)\big(X_{t-h,j} - \bar X_j\big),
$$

using the $n - |h|$ overlapping pairs and a divisor of $n-|h|$. Each entry is
then pushed through the inverse link,

$$
\hat\gamma^{Z}_{ij}(h) = L_{ij}^{-1}\big(\hat\gamma^{X}_{ij}(h)\big).
$$

The inverse is computed numerically: $L_{ij}$ is evaluated on a grid of $u$
values in $[-1, 1]$ (fine near $\pm 1$, where the link is steepest, coarse in
the middle), a natural cubic spline is fitted through the pairs
$(L_{ij}(u_g), u_g)$, and the spline is evaluated at the observed covariance.
Values outside the range of the link are clamped to $\pm 1$ (with a warning).

Collecting the entries gives the two matrices that everything downstream uses:

$$
S_0 = \big[\hat\gamma^Z_{ij}(0)\big]_{ij} \quad (\text{lag-0 latent covariance}), \qquad
S_{10} = \big[\hat\gamma^Z_{ij}(1)\big]_{ij} = \widehat{\operatorname{Cov}}(Z_t, Z_{t-1}) \quad (\text{lag-1, current-by-past}).
$$

The lag $-1$ matrix is the transpose, $\widehat{\operatorname{Cov}}(Z_{t-1}, Z_t) = S_{10}^\top$.

### 2.5 Positive-definiteness of the joint latent covariance

The entrywise inverse link does not preserve positive semidefiniteness: even
though the observed lag-0 covariance is PSD, the matrix of latent estimates need
not be, and more importantly the *joint* covariance of $(Z_{t-1}, Z_t)$,

$$
J = \begin{pmatrix} S_0 & S_{10}^\top \\ S_{10} & S_0 \end{pmatrix},
$$

need not be. This matters because (Section 4) the innovation covariance implied
by any transition matrix $B$ is a congruence of $J$,

$$
S_\varepsilon(B) = \begin{pmatrix} -B & I \end{pmatrix} J \begin{pmatrix} -B & I \end{pmatrix}^\top,
$$

so $J \succeq 0$ guarantees $S_\varepsilon(B) \succeq 0$ for *every* $B$, while
an indefinite $J$ can make the Gaussian likelihood unbounded. With
`pd_approx = TRUE` the constructor repairs $J$ once, at construction, by
alternating projections between the PSD cone (via `Matrix::nearPD`) and the
set of matrices with the block-Toeplitz structure above (equal diagonal blocks,
transpose-consistent off-diagonal blocks), both of which are convex sets. A
fallback shrinks $S_{10}$ toward zero, whose block-diagonal limit is PD once
$S_0$ is. The repaired $S_0, S_{10}$ replace the originals. With the default
`pd_approx = FALSE`, the multitask fit checks each subject's $J$ and warns if it
is not PSD.

### 2.6 The single-subject estimator, for orientation

If one only wanted $B$, the Yule–Walker equation $\Gamma(1) = B\,\Gamma(0)$
(which follows from multiplying the VAR equation by $Z_{t-1}^\top$ and taking
expectations) gives the moment estimator

$$
\hat B_{\text{YW}} = S_{10}\, S_0^{-1}.
$$

This is what `fit_timecop()` returns. Everything in this document generalizes
this idea to $K$ subjects with shared structure and sparsity penalties, and
adds the precision network.

---

## 3. The multitask graphical VAR model

Now return to $K$ subjects. Each subject $k$ has its own latent VAR(1),

$$
Z^{(k)}_t = B_k\, Z^{(k)}_{t-1} + \varepsilon^{(k)}_t, \qquad \varepsilon^{(k)}_t \sim \mathcal{N}_d(0, \Omega_k^{-1}),
$$

with its own observation transforms as in Section 2.2 (same families across
subjects; marginal parameters estimated per subject).

### 3.1 Common-plus-unique decompositions

The multitask model assumes each subject's two networks are perturbations of a
shared network:

$$
\boxed{\;B_k = \mu + \Delta_k\;} \qquad\qquad
\boxed{\;\Omega_k = M + E_k, \qquad \operatorname{diag}(M) = 0,\quad \Omega_k \succ 0\;}
$$

* $\mu \in \mathbb{R}^{d\times d}$ is the **common temporal network** and $\Delta_k$ is subject $k$'s **unique temporal deviation**. All $d^2$ entries of each are free.
* $M \in \mathbb{R}^{d\times d}$ is the **common contemporaneous network**, symmetric with zero diagonal, so it carries only shared *edges*. $E_k$ is symmetric and carries subject $k$'s **unique edges** *and* its entire precision diagonal. Fixing $\operatorname{diag}(M) = 0$ resolves the diagonal ambiguity of the split and leaves each subject's innovation variances unpenalized and unconstrained.

The word "multitask" refers to the structure of multi-task learning: a shared
parameter that borrows strength across tasks, plus task-specific parameters
that are shrunk toward zero.

### 3.2 Identifiability

Without further constraints, $(\mu, \{\Delta_k\})$ is not identified from
$\{B_k\}$: for any $D$, $(\mu + D, \{\Delta_k - D\})$ gives the same $B_k$. The
same holds for $(M, \{E_k\})$ on the off-diagonal. The decomposition becomes
meaningful only through the sparsity penalties (Section 5), which favor
representations in which the shared part absorbs what is common across subjects
and the deviations are sparse. In particular:

* if $\lambda_\Delta = 0$ the deviations absorb everything and any $\lambda_\mu > 0$ drives $\mu$ to zero, so $\mu$ only carries signal when deviations are penalized;
* the estimated $\mu$ should be read as "the sparse part shared by all subjects", not as the cross-subject mean of the $B_k$.

### 3.3 Innovation covariance implied by a transition matrix

The two blocks of parameters are coupled through the innovation covariance.
Given subject $k$'s latent moments, define the model-implied innovation
covariance at a transition matrix $B$ as

$$
S_{\varepsilon,k}(B) = C_k - B H_k^\top - H_k B^\top + B\, G_k\, B^\top,
$$

where $C_k$, $H_k$ and $G_k$ are the (population or sample) second moments of
the current vector, the current-by-past pair, and the past vector,
respectively. Under stationarity $C_k = G_k = \Gamma_k(0)$ and
$H_k = \Gamma_k(1)$, so the package uses

$$
C_k = G_k = S_{0,k}, \qquad H_k = S_{10,k}.
$$

At the true parameters, $\Gamma(1) = B\,\Gamma(0)$ and
$\Gamma(0) = B\,\Gamma(0)B^\top + \Sigma$, so
$S_{\varepsilon}(B) = \Gamma(0) - B\Gamma(0)B^\top = \Sigma$, as it should. At
the sample level, $S_{\varepsilon,k}(B)$ is exactly the residual covariance
$\tfrac{1}{N_k}\sum_t (z_t - B z_{t-1})(z_t - B z_{t-1})^\top$ one would obtain
by regressing the latent series on its lag, had the latent series been
observed. This is the sense in which the whole estimator is *covariance-only*:
it never needs the latent series, only $S_{0,k}$ and $S_{10,k}$.

---

## 4. From data to the covariance-only objective

### 4.1 Gaussian conditional log-likelihood of a VAR(1)

Suppose for a moment that the latent series $z^{(k)}_1, \dots, z^{(k)}_{n_k}$
were observed. Conditional on $z^{(k)}_1$, the Gaussian log-likelihood of
subject $k$ is

$$
\ell_k(B_k, \Omega_k) = \frac{N_k}{2}\log\det\Omega_k - \frac{1}{2}\sum_{t=2}^{n_k} \big(z_t - B_k z_{t-1}\big)^\top \Omega_k \big(z_t - B_k z_{t-1}\big) + \text{const}.
$$

Writing the quadratic form as a trace and pulling out the residual covariance,

$$
\sum_{t} (z_t - B_k z_{t-1})^\top \Omega_k (z_t - B_k z_{t-1}) = N_k \operatorname{tr}\!\big(S_{\varepsilon,k}(B_k)\,\Omega_k\big),
$$

so the negative log-likelihood is

$$
-\ell_k(B_k, \Omega_k) = \frac{N_k}{2}\Big[-\log\det\Omega_k + \operatorname{tr}\!\big(S_{\varepsilon,k}(B_k)\,\Omega_k\big)\Big] + \text{const}.
$$

This depends on the data only through the second moments $C_k, H_k, G_k$
inside $S_{\varepsilon,k}$.

### 4.2 The plug-in step

Since the latent series is not observed, the package replaces the unavailable
latent sample moments by the copula estimates of Section 2.4:

$$
C_k \leftarrow S_{0,k}, \qquad G_k \leftarrow S_{0,k}, \qquad H_k \leftarrow S_{10,k}.
$$

The resulting function of $(B_k, \Omega_k)$ is a Gaussian *pseudo*-likelihood:
it has the algebraic form of the latent Gaussian likelihood, evaluated at
method-of-moments estimates of the latent second moments. This is the same
logic that turns the population Yule–Walker equation into the single-subject
estimator, extended to a full likelihood so that a precision matrix can be
estimated and penalties can be attached.

Two consequences of the plug-in are worth keeping in mind:

* Subjects are weighted by $N_k$: a longer series contributes proportionally more to every objective below, exactly as in pooled least squares.
* The latent scale is a convention (Section 2.1), and the inverse link returns latent quantities on that unit-variance scale at both lags, so $S_{0,k}$ and $S_{10,k}$ are mutually consistent. With discrete marginals the estimated lag-0 diagonal is not pinned to one (for Bernoulli it equals one algebraically; for Poisson it is close to one and shares its sampling noise with the lag-1 entries, so the Yule–Walker ratio largely cancels it).

---

## 5. The penalized objective

Collecting the $K$ pseudo-likelihoods and attaching one lasso-type penalty to
each of the four networks, the multitask graphical VAR estimator is the
minimizer of

$$
\boxed{
\begin{aligned}
F(\mu, \{\Delta_k\}, M, \{E_k\}) \;=\; & \sum_{k=1}^K \frac{N_k}{2}\Big[-\log\det\Omega_k + \operatorname{tr}\!\big(S_{\varepsilon,k}(\mu + \Delta_k)\,\Omega_k\big)\Big] \\
& + \lambda_\mu \|W_\mu \circ \mu\|_1 + \lambda_\Delta \sum_{k=1}^K \|W_{\Delta,k} \circ \Delta_k\|_1 \\
& + \lambda_M \|W_M \circ M\|_{1,\mathrm{off}} + \lambda_E \sum_{k=1}^K \|W_{E,k} \circ E_k\|_{1,\mathrm{off}}
\end{aligned}}
$$

subject to $\Omega_k = M + E_k$, $\operatorname{diag}(M) = 0$, $M = M^\top$,
$E_k = E_k^\top$, and $\Omega_k \succ 0$. Here $\circ$ is the entrywise
(Hadamard) product and the $W$'s are nonnegative weight matrices: all ones for
the plain lasso, or data-driven adaptive weights (Section 10). Four tuning
parameters control sparsity:

| Penalty | Acts on | Effect when large |
|---|---|---|
| $\lambda_\mu$ | common temporal network $\mu$ | fewer shared lagged effects |
| $\lambda_\Delta$ | unique temporal deviations $\Delta_k$ | subjects' dynamics pulled toward $\mu$ |
| $\lambda_M$ | common contemporaneous edges $M$ | fewer shared partial correlations |
| $\lambda_E$ | unique contemporaneous edges $E_k$ | subjects' networks pulled toward $M$ |

The temporal penalties act on all $d^2$ entries (including the autoregressive
diagonal), while the precision penalties act on off-diagonal entries only, so
innovation variances are never shrunk.

### 5.1 Structure of the objective

Write the objective as $F = f(\mu, \Delta, \Omega) + P_{\text{temp}}(\mu, \Delta) + P_{\text{prec}}(M, E)$, where $f$ is the smooth pseudo-likelihood term. Two facts drive the algorithm:

1. **For fixed $\{\Omega_k\}$, $F$ is a convex quadratic plus an $\ell_1$ penalty in $(\mu, \{\Delta_k\})$.** Expanding the trace and dropping terms that do not involve $B_k$,
   $$
   \frac{N_k}{2}\operatorname{tr}\!\big(S_{\varepsilon,k}(B_k)\Omega_k\big) = \frac{N_k}{2}\operatorname{tr}\!\big(B_k S_{0,k} B_k^\top \Omega_k\big) - N_k \operatorname{tr}\!\big(B_k S_{10,k}^\top \Omega_k\big) + \text{const},
   $$
   using $\operatorname{tr}(B H^\top \Omega) = \operatorname{tr}(H B^\top \Omega)$ for symmetric $\Omega$. This is a *precision-weighted* (generalized) least-squares problem, solved in Section 7.

2. **For fixed $\{B_k\}$, $F$ is a multitask graphical lasso in $(M, \{E_k\})$.** The residual covariances $S_{\varepsilon,k}(B_k)$ are fixed matrices, and the objective in $\Omega_k$ is the weighted Gaussian graphical-model negative log-likelihood $w_k[-\log\det\Omega_k + \operatorname{tr}(S_{\varepsilon,k}\Omega_k)]$ with $w_k = N_k/2$, coupled across subjects by the shared $M$. This is convex in $\{\Omega_k\}$ and solved in Section 8.

Each block is convex, but $F$ is not jointly convex in $(B, \Omega)$ (the
coupling $\operatorname{tr}(S_\varepsilon(B)\,\Omega)$ is bilinear in $B B^\top$
and $\Omega$). The natural algorithm is therefore block coordinate descent,
which is what the package does.

---

## 6. Optimization I: the outer alternation

At fixed $(\lambda_\mu, \lambda_\Delta, \lambda_M, \lambda_E)$ the package
minimizes $F$ by alternating exact (to tolerance) minimizations over the two
blocks:

> **Algorithm A (alternation at fixed penalties)**
>
> *Input:* $S_{0,k}$, $S_{10,k}$, $N_k$ for $k = 1..K$; penalties; weights; optional warm starts.
>
> 0. **Initialize the temporal block.** If warm starts $(\mu^{(0)}, \Delta^{(0)})$ are supplied, use them. Otherwise run the temporal solver with $\Omega_k = I$ (an unweighted, i.e. ordinary least-squares, multitask fit) to obtain $(\mu^{(0)}, \Delta^{(0)})$.
> 1. **For** $s = 1, 2, \dots$ until convergence:
>    1. Form $B_k = \mu^{(s-1)} + \Delta_k^{(s-1)}$ and $S_{\varepsilon,k} = S_{\varepsilon,k}(B_k)$.
>    2. **Precision block.** Run the multitask graphical-lasso ADMM (Section 8) on $\{S_{\varepsilon,k}\}$ with weights $w_k = N_k/2$, warm-started from the previous ADMM state $(\Omega, M, E, U)$, to obtain $\Omega_k^{(s)}, M^{(s)}, E_k^{(s)}$.
>    3. **Temporal block.** Run the precision-weighted FISTA solver (Section 7) with $\Omega_k = \Omega_k^{(s)}$, warm-started from $(\mu^{(s-1)}, \Delta^{(s-1)})$, to obtain $\mu^{(s)}, \Delta_k^{(s)}$.
>    4. Stop if $\max\big(\|\mu^{(s)} - \mu^{(s-1)}\|_\infty,\ \max_k \|\Delta_k^{(s)} - \Delta_k^{(s-1)}\|_\infty\big) < \text{tol}$ (default $10^{-5}$), or after `max_outer` (default 50) iterations.
> 2. Return $\mu, \Delta_k, B_k, \Omega_k, M, E_k$, the ADMM duals $U_k$ (for warm-starting the next call), and the final value of $F$.

Because each block is solved to (near) optimality and both subproblems are
convex, the sequence of objective values is monotonically non-increasing, so
the alternation converges to a stationary point of $F$. The outer tolerance is
deliberately looser than the inner solver tolerances: using the FISTA
tolerance ($10^{-7}$) as the outer tolerance makes the alternation run to its
cap for no gain.

Note the order within an outer iteration: the precision block goes first,
because on entry the temporal parameters are the only thing available, and the
precision block is what turns them into a first estimate of $\Omega_k$.

If $\lambda_M$ or $\lambda_E$ is not supplied to the alternation, it defaults
to one tenth of the corresponding KKT anchor (Section 9) computed at the
initial $B_k$.

---

## 7. Optimization II: the temporal block (FISTA)

### 7.1 The subproblem

Fix $\{\Omega_k\}$. Using the expansion of Section 5.1, the temporal
subproblem is

$$
\min_{\mu, \{\Delta_k\}} \; \sum_{k=1}^K \Big[\tfrac{1}{2}\operatorname{tr}\!\big(B_k\, \mathsf{G}_k\, B_k^\top \Omega_k\big) - \operatorname{tr}\!\big(B_k\, \mathsf{H}_k^\top \Omega_k\big)\Big] + \lambda_\mu \|W_\mu \circ \mu\|_1 + \lambda_\Delta \sum_k \|W_{\Delta,k} \circ \Delta_k\|_1,
$$

with $B_k = \mu + \Delta_k$ and the $N_k$-weighted moments

$$
\mathsf{G}_k = N_k\, S_{0,k}, \qquad \mathsf{H}_k = N_k\, S_{10,k}.
$$

The solver is agnostic to this scaling; the wrapper `fit_multitask()` forms
$\mathsf{G}_k$ and $\mathsf{H}_k$ so that the temporal objective is on the same
scale as the pseudo-likelihood and the EBIC. The objective is a smooth convex
function $f$ (a quadratic in the stacked parameters) plus a separable
non-smooth penalty $g$, the textbook setting for proximal gradient methods.

### 7.2 Gradient

Using $\partial \operatorname{tr}(B \mathsf{G} B^\top \Omega)/\partial B = 2\,\Omega B \mathsf{G}$ and $\partial \operatorname{tr}(B \mathsf{H}^\top \Omega)/\partial B = \Omega \mathsf{H}$ for symmetric $\mathsf{G}, \Omega$, the gradient of the smooth part with respect to each subject's total matrix is

$$
\nabla_{B_k} f = \Omega_k \big(B_k \mathsf{G}_k - \mathsf{H}_k\big) = N_k\, \Omega_k \big(B_k S_{0,k} - S_{10,k}\big).
$$

The term in parentheses is the Yule–Walker residual for subject $k$; it
vanishes at $\hat B_{\text{YW},k} = S_{10,k} S_{0,k}^{-1}$ regardless of
$\Omega_k$, which is why at zero penalty the multitask solver reproduces the
per-subject Yule–Walker estimates exactly (Section 13). The precision matrix
enters only as a left-multiplication: it reweights the residual directions, so
that a residual in a low-variance innovation direction costs more. This is the
familiar fact that generalized least squares differs from ordinary least
squares only through the weighting of the estimating equations.

By the chain rule for $B_k = \mu + \Delta_k$,

$$
\nabla_{\Delta_k} f = \nabla_{B_k} f =: \mathsf{g}_k, \qquad \nabla_{\mu} f = \sum_{k=1}^K \mathsf{g}_k .
$$

### 7.3 Stacked representation and the Lipschitz constant

Stack the $K + 1$ parameter blocks into one tall matrix
$\Theta = [\mu; \Delta_1; \dots; \Delta_K] \in \mathbb{R}^{(K+1)d \times d}$ (the
implementation stores the transposes, see Section 15; this is immaterial for an
entrywise penalty). The smooth part is then $f(\Theta) = \tfrac12 \langle \Theta, \mathcal{Q}\Theta\rangle - \langle \Theta, R\rangle$ for a symmetric positive semidefinite linear operator $\mathcal{Q}$ on the stacked space, whose action is

$$
(\mathcal{Q}\Theta)_k = \Omega_k(\mu + \Delta_k)\mathsf{G}_k \quad (k \ge 1), \qquad (\mathcal{Q}\Theta)_0 = \sum_k (\mathcal{Q}\Theta)_k .
$$

The gradient is Lipschitz with constant $L = \lambda_{\max}(\mathcal{Q})$. The
operator is never formed as a matrix (it would be $(K+1)d^2 \times (K+1)d^2$);
instead $L$ is estimated by **power iteration** on the map
$\Theta \mapsto \mathcal{Q}\Theta$: starting from a random unit-norm $\Theta$,
repeat $\Theta \leftarrow \mathcal{Q}\Theta / \|\mathcal{Q}\Theta\|_F$ and take
$L = \langle \Theta, \mathcal{Q}\Theta\rangle$ at convergence (relative change
below $10^{-6}$, at most 100 iterations). This is done once per call.

### 7.4 Proximal operator

The non-smooth part is $g(\Theta) = \sum_{\text{entries}} \tau_{ab}\lvert\Theta_{ab}\rvert$ with a block-structured threshold matrix

$$
\tau = \begin{bmatrix} \lambda_\mu W_\mu \\ \lambda_\Delta W_{\Delta,1} \\ \vdots \\ \lambda_\Delta W_{\Delta,K}\end{bmatrix},
$$

whose proximal operator with step $1/L$ is entrywise soft-thresholding,
$\operatorname{prox}_{g/L}(V) = \mathcal{S}_{\tau/L}(V)$. A plain proximal
gradient (ISTA) step is therefore

$$
\Theta^+ = \mathcal{S}_{\tau/L}\Big(\Theta - \tfrac{1}{L}\nabla f(\Theta)\Big).
$$

### 7.5 FISTA with monotone restart

The solver uses Nesterov acceleration in the form of FISTA (Beck and Teboulle,
2009) with a monotonicity safeguard:

> **Algorithm B (temporal solver)**
>
> Initialize $\Theta^{(0)}$ (zeros, or the warm start), $Y^{(0)} = \Theta^{(0)}$, momentum scalar $q_0 = 1$.
> For $i = 1, 2, \dots$:
> 1. Candidate: $\tilde\Theta = \mathcal{S}_{\tau/L}\big(Y^{(i-1)} - \nabla f(Y^{(i-1)})/L\big)$.
> 2. **If** $F_{\text{temp}}(\tilde\Theta) > F_{\text{temp}}(\Theta^{(i-1)})$ (the accelerated step went uphill): discard it, take an ISTA step from the current iterate, $\Theta^{(i)} = \mathcal{S}_{\tau/L}\big(\Theta^{(i-1)} - \nabla f(\Theta^{(i-1)})/L\big)$, and reset the momentum: $q_i = 1$, $Y^{(i)} = \Theta^{(i)}$.
>    **Else** accept $\Theta^{(i)} = \tilde\Theta$ and update the momentum: $q_i = \big(1 + \sqrt{1 + 4 q_{i-1}^2}\big)/2$, $Y^{(i)} = \Theta^{(i)} + \frac{q_{i-1} - 1}{q_i}\big(\Theta^{(i)} - \Theta^{(i-1)}\big)$.
> 3. Stop when the relative change in the penalized objective is below `tol` (default $10^{-7}$) or after `max_iter` (default 1000) iterations.

FISTA attains the optimal $O(1/i^2)$ rate for this problem class; the restart
(in the spirit of O'Donoghue and Candès, 2015) prevents the oscillations that
plain FISTA exhibits near the solution and makes the objective sequence
non-increasing, which matters because the outer alternation relies on
monotone decrease. Warm starts are essential in the outer loop and along the
penalty path: a solve from a nearby solution typically converges in a handful
of iterations.

---

## 8. Optimization III: the precision block (ADMM)

### 8.1 The subproblem

Fix $\{B_k\}$ and set $S_k := S_{\varepsilon,k}(B_k)$, $w_k := N_k/2$. The
precision subproblem is

$$
\min_{\{\Omega_k\}, M, \{E_k\}} \; \sum_{k=1}^K w_k\Big[-\log\det\Omega_k + \operatorname{tr}(S_k\Omega_k)\Big] + \lambda_M\|W_M \circ M\|_{1,\mathrm{off}} + \lambda_E\sum_k \|W_{E,k} \circ E_k\|_{1,\mathrm{off}}
$$

$$
\text{subject to}\quad \Omega_k = M + E_k \ \ (k = 1..K), \qquad \operatorname{diag}(M) = 0, \qquad \Omega_k \succ 0 .
$$

This is a *joint* graphical lasso in the sense of Danaher, Wang and Witten
(2014), but with a common-plus-unique parameterization of the precision
matrices rather than a fused or group penalty. The equality constraints are
exactly the form ADMM (Boyd et al., 2011) is designed for: one set of variables
($\Omega_k$) carries the smooth likelihood and the positive-definiteness
constraint, the other set ($M, E_k$) carries the penalties, and they are tied
by a linear constraint.

### 8.2 Scaled augmented Lagrangian

Introduce scaled dual variables $U_k \in \mathbb{R}^{d\times d}$ and a penalty
parameter $\rho > 0$ (default 1). The scaled-form augmented Lagrangian is

$$
\mathcal{L}_\rho = \sum_k w_k\big[-\log\det\Omega_k + \operatorname{tr}(S_k\Omega_k)\big] + P_{\text{prec}}(M, E) + \frac{\rho}{2}\sum_k \big\|\Omega_k - M - E_k + U_k\big\|_F^2 + \text{const}.
$$

ADMM cycles through three steps: minimize over $\{\Omega_k\}$, minimize over
$(M, \{E_k\})$, then update the duals. Each step has a closed form or a cheap
inner solver.

### 8.3 Step 1: per-subject precision update (closed form)

For each $k$, with $A_k := M + E_k - U_k$, solve

$$
\min_{\Omega \succ 0}\; w_k\big[-\log\det\Omega + \operatorname{tr}(S_k\Omega)\big] + \frac{\rho}{2}\|\Omega - A_k\|_F^2 .
$$

The stationarity condition is $w_k(S_k - \Omega^{-1}) + \rho(\Omega - A_k) = 0$, i.e.

$$
\Omega - \frac{w_k}{\rho}\,\Omega^{-1} = A_k - \frac{w_k}{\rho}\,S_k =: C_k .
$$

Symmetrize $C_k$ and eigendecompose it, $C_k = Q\,\operatorname{diag}(c_1, \dots, c_d)\,Q^\top$. Because the left-hand side is a matrix function of $\Omega$, the solution shares the eigenvectors $Q$, and each eigenvalue $\theta_j$ solves the scalar quadratic $\theta_j - (w_k/\rho)/\theta_j = c_j$, whose positive root is

$$
\theta_j = \frac{c_j + \sqrt{c_j^2 + 4 w_k/\rho}}{2} > 0 .
$$

So $\Omega_k \leftarrow Q\,\operatorname{diag}(\theta)\,Q^\top$ is automatically
symmetric positive definite, whatever $A_k$ and $S_k$ are; the PD constraint
never has to be enforced separately. This is the same update as in the ADMM
for the single-task graphical lasso, with the likelihood weight $w_k$ carried
through.

### 8.4 Step 2: common-plus-unique decomposition (edgewise)

With $R_k := \Omega_k + U_k$ fixed, solve

$$
\min_{M, \{E_k\}}\; \lambda_M\|W_M \circ M\|_{1,\mathrm{off}} + \lambda_E\sum_k\|W_{E,k} \circ E_k\|_{1,\mathrm{off}} + \frac{\rho}{2}\sum_k\|R_k - M - E_k\|_F^2, \qquad \operatorname{diag}(M) = 0 .
$$

This problem separates completely across matrix positions.

* **Diagonal.** $M_{ii} = 0$ by constraint and there is no penalty, so $(E_k)_{ii} = (R_k)_{ii}$.
* **Off-diagonal pair $(i, j)$, $i < j$.** Writing $a_k = (R_k)_{ij}$, $m = M_{ij}$, $e_k = (E_k)_{ij}$ (symmetry means the $(j,i)$ entries are the same unknowns and contribute an identical term, which doubles both the loss and the penalty and therefore cancels), the problem is the scalar common-plus-unique lasso
  $$
  \min_{m, e_1..e_K}\; \lambda_M W_{M,ij}\lvert m\rvert + \lambda_E\sum_k W_{E,k,ij}\lvert e_k\rvert + \frac{\rho}{2}\sum_k (a_k - m - e_k)^2 .
  $$
  It is convex with a separable non-smooth part, so coordinate descent converges to its global minimizer. Both coordinate updates are soft-thresholds:
  $$
  e_k \leftarrow \mathcal{S}_{\lambda_E W_{E,k,ij}/\rho}\big(a_k - m\big), \qquad
  m \leftarrow \mathcal{S}_{\lambda_M W_{M,ij}/(\rho K)}\Big(\frac{1}{K}\sum_k (a_k - e_k)\Big).
  $$
  The $m$-update follows from $\frac{\rho}{2}\sum_k (a_k - e_k - m)^2 = \frac{\rho K}{2}\big(m - \overline{a - e}\big)^2 + \text{const}$. The two updates are iterated (at most 100 times, to a tolerance of $10^{-8}$) from the previous ADMM iterate's $(m, e)$.

The interpretation is transparent: a shared edge appears in $M$ only if the
*average* signal across subjects survives a threshold that shrinks like $1/K$
(evidence accumulates across subjects), while a unique edge appears in $E_k$
only if subject $k$'s residual signal *after removing the shared part* survives
its own threshold.

### 8.5 Step 3: dual update and stopping

$$
U_k \leftarrow U_k + \Omega_k - M - E_k, \qquad k = 1..K .
$$

Convergence is declared using the standard primal/dual residual test of Boyd
et al. (2011, Section 3.3). With primal residual
$r = \big(\sum_k\|\Omega_k - M - E_k\|_F^2\big)^{1/2}$ and dual residual
$s = \rho\big(\sum_k\|(M + E_k) - (M^{\text{old}} + E_k^{\text{old}})\|_F^2\big)^{1/2}$, stop when

$$
r \le \sqrt{Kd^2}\,\epsilon_{\text{abs}} + \epsilon_{\text{rel}}\max\Big(\big(\textstyle\sum_k\|\Omega_k\|_F^2\big)^{1/2}, \big(\sum_k\|M + E_k\|_F^2\big)^{1/2}\Big), \qquad
s \le \sqrt{Kd^2}\,\epsilon_{\text{abs}} + \epsilon_{\text{rel}}\,\rho\,\big(\textstyle\sum_k\|U_k\|_F^2\big)^{1/2},
$$

with defaults $\epsilon_{\text{abs}} = 10^{-5}$, $\epsilon_{\text{rel}} = 10^{-4}$, and at most 1000 iterations.

### 8.6 Initialization and warm starts

Cold start: $M = 0$, $E_k = \operatorname{diag}(1/\operatorname{diag}(S_k))$
(the precision of independent innovations), $\Omega_k = E_k$, $U_k = 0$. Inside
the outer alternation the ADMM is warm-started from the previous outer
iteration's full state $(\Omega, M, E, U)$, and across penalty-grid points from
the previous grid point's state. Because consecutive subproblems differ only
slightly, a warm-started ADMM typically converges in one or a few iterations;
this is where most of the speed of the grid search comes from.

### 8.7 What is returned

ADMM returns two closely related objects: the $\Omega_k$ from Step 1, which are
exactly positive definite but only approximately equal to $M + E_k$ (the gap is
the primal residual, of order $10^{-5}$), and the sparse decomposition
$(M, E_k)$, whose entries are exact zeros after soft-thresholding. The package
uses $\Omega_k$ wherever a likelihood is evaluated (the EBIC, the
precision-weighted temporal gradient) and $(M, E_k)$ wherever sparsity is read
off (degrees of freedom, edge sets). Edge sets should therefore be read from
`M_hat` and `E_hat`, not from `Omega_hat`, whose "zeros" are only
residual-sized.

---

## 9. Penalty grids and KKT anchors

The four penalties are selected by grid search (Section 12), so the grid must
bracket the range in which each penalty actually matters. For the lasso this
range has a known top: the smallest penalty at which the fully sparse solution
is optimal, the analogue of $\lambda_{\max} = \max\lvert X^\top y\rvert$ in
lasso regression. The package derives this anchor separately for each of the
four axes from the Karush–Kuhn–Tucker (KKT) conditions, then descends three
decades on a log scale:

$$
\lambda^{(1)} = \lambda_{\max}, \quad \lambda^{(n)} = 10^{-3}\lambda_{\max}, \quad \text{log-equispaced in between.}
$$

A user-supplied vector bypasses the automatic grid for that axis; a supplied
scalar collapses the axis to a single value.

### 9.1 Temporal anchors

For a convex problem $\min_\theta f(\theta) + \lambda\|W \circ \theta\|_1$, the
zero vector is optimal if and only if $\lvert\nabla f(0)_{ab}\rvert \le \lambda W_{ab}$ for all entries. In the temporal subproblem with $\Omega_k = I$ (the state at initialization), the gradient at $\mu = \Delta_k = 0$ is $\nabla_{\Delta_k} f = -\mathsf{H}_k = -N_k S_{10,k}$ and $\nabla_\mu f = -\sum_k \mathsf{H}_k$. Hence

$$
\lambda_{\mu,\max} = \max_{ij}\frac{\big\lvert\sum_k N_k (S_{10,k})_{ij}\big\rvert}{(W_\mu)_{ij}}, \qquad
\lambda_{\Delta,\max} = \max_k\max_{ij}\frac{N_k\lvert(S_{10,k})_{ij}\rvert}{(W_{\Delta,k})_{ij}} .
$$

(With uniform weights the denominators are one.) These carry the factor $N_k$
because the loss does.

### 9.2 Precision anchors

In the precision subproblem, consider the fully sparse candidate: diagonal
$\Omega_k$, i.e. all off-diagonals of $M$ and $E_k$ zero. The gradient of the
smooth part with respect to an off-diagonal entry of $\Omega_k$ is
$w_k(S_k - \Omega_k^{-1})_{ij}$, and at a diagonal $\Omega_k$ the inverse is
diagonal, so the gradient is $w_k (S_k)_{ij}$. The KKT condition for keeping
$(E_k)_{ij}$ at zero is $\lvert w_k (S_k)_{ij}\rvert \le \lambda_E W_{E,k,ij}$;
for the shared entry $M_{ij}$, whose gradient is the sum over subjects, it is
$\lvert\sum_k w_k (S_k)_{ij}\rvert \le \lambda_M W_{M,ij}$. Hence

$$
\lambda_{E,\max} = \max_k\max_{i\ne j}\frac{w_k\lvert(S_{\varepsilon,k})_{ij}\rvert}{W_{E,k,ij}}, \qquad
\lambda_{M,\max} = \max_{i\ne j}\frac{\big\lvert\sum_k w_k (S_{\varepsilon,k})_{ij}\big\rvert}{W_{M,ij}}, \qquad w_k = \frac{N_k}{2}.
$$

The residual covariances used here are the ones implied by the per-subject
Yule–Walker transitions, $S_{\varepsilon,k}(\hat B_{\text{YW},k})$: the anchor
shifts slightly as $B_k$ is updated during the alternation, which is
immaterial on a log-spaced grid.

Two features of these anchors are easy to get wrong and are worth stating
explicitly (both were bugs in earlier versions of the code and are documented
in `multitask_precision_grid_diagnosis.md`):

* The anchors are proportional to $w_k = N_k/2$ because the penalty must compete with an $N_k/2$-weighted likelihood. Anchoring on the unweighted covariance, as one would for a single-subject graphical lasso with an unweighted likelihood, puts the whole grid two orders of magnitude too low and makes every candidate model dense.
* Under adaptive weights the gradients must be divided entrywise by the weights before taking the maximum. Precision weights are typically larger than one (pilot precision entries are below one in magnitude), so the weighted problem's active range sits one to three decades below the unweighted anchor and is stretched; without the correction the EBIC selects at the grid boundary.

For the same reason the temporal division uses the weights in the same
orientation as the gradient; the implementation stores the transposed
parameterization and therefore divides by $W^\top$ (Section 15).

### 9.3 Boundary diagnostics

Whenever the selected value of an automatically constructed axis falls at
either end of its grid, the fit emits a message. A selection at the top means
"fully sparse was best", which is legitimate when that block is truly empty; a
selection at the bottom is the classic sign of a grid that does not reach far
enough.

---

## 10. Adaptive weights

The lasso shrinks large and small coefficients alike, which biases the large
ones and tends to over-select. The adaptive lasso (Zou, 2006) fixes this with
entry-specific weights computed from a pilot estimate,

$$
W_{ab} = \frac{1}{(\lvert\hat\theta^{\text{pilot}}_{ab}\rvert + \epsilon)^{\gamma}}, \qquad \epsilon = 10^{-3}\max_{ab}\lvert\hat\theta^{\text{pilot}}_{ab}\rvert,
$$

with $\gamma = 1$ by default (`penalty_gamma`). Entries that are large in the
pilot get small weights and are barely penalized; entries near zero get large
but finite weights. The $\epsilon$ term keeps weights finite so that an exact
zero in a penalized pilot does not permanently lock an edge out.

With `penalty = "adaptive"`, all four weight sets are built from one
per-subject pilot $(\hat B_k, \hat\Omega_k)$, split the same way the estimator
splits its parameters:

$$
\hat\mu^{0} = \frac{1}{K}\sum_k \hat B_k, \qquad \hat\Delta^{0}_k = \hat B_k - \hat\mu^{0}, \qquad
\hat M^{0} = \operatorname{offdiag}\Big(\frac{1}{K}\sum_k \hat\Omega_k\Big), \qquad \hat E^{0}_k = \hat\Omega_k - \hat M^{0},
$$

and $W_\mu, W_{\Delta,k}, W_M, W_{E,k}$ are the reciprocal-magnitude weights of
these four quantities (precision weights on off-diagonals only, with zero
diagonal). The pilot itself is chosen by `adaptive_weights_from`:

* `"unpenalized"` (default): $\hat B_k = S_{10,k}S_{0,k}^{-1}$ (Yule–Walker) and $\hat\Omega_k = \big(S_{\varepsilon,k}(\hat B_k) + 10^{-3}I\big)^{-1}$. Because GLS and OLS coincide in an unpenalized VAR with the same regressors in every equation, this is the unpenalized single-subject graphical VAR.
* `"penalized"`: a per-subject EBIC-selected single-task graphical VAR fit (`fit_graphical_var()`) on a small grid.

The weights are computed before the penalty grids so that the grids can be
calibrated to the weighted KKT conditions (Section 9). With
`penalty = "lasso"` all weights are one.

`penalty = "scad"` (Fan and Li, 2001), implemented via local linear
approximation for the ordinary multitask VAR, is not implemented for the
graphical-VAR path: with `gvar = TRUE` it currently falls back to the plain
lasso.

---

## 11. Model selection: the joint EBIC

The four penalties are chosen to minimize the extended Bayesian information
criterion (Chen and Chen, 2008; Foygel and Drton, 2010) of the complete model.
For a fitted model with $\hat\mu, \hat\Delta_k, \hat B_k, \hat\Omega_k, \hat M, \hat E_k$:

**Log-likelihood.** The Gaussian pseudo-log-likelihood of Section 4, evaluated at the estimates and dropping constants,

$$
\hat\ell = \sum_{k=1}^K \frac{N_k}{2}\Big[\log\det\hat\Omega_k - \operatorname{tr}\!\big(S_{\varepsilon,k}(\hat B_k)\,\hat\Omega_k\big)\Big].
$$

**Degrees of freedom.** Nonzero counts (entries with magnitude above $10^{-8}$) in each network:

$$
\text{df}_\mu = \#\{\hat\mu_{ij} \ne 0\}, \quad
\text{df}_\Delta = \sum_k \#\{(\hat\Delta_k)_{ij} \ne 0\}, \quad
\text{df}_M = \#\{\hat M_{ij} \ne 0,\ i < j\}, \quad
\text{df}_E = \sum_k \#\{(\hat E_k)_{ij} \ne 0,\ i < j\},
$$

and $\text{df} = \text{df}_\mu + \text{df}_\Delta + \text{df}_M + \text{df}_E$.
Temporal networks count all entries (the autoregressive diagonal is a
parameter); precision networks count upper-triangular off-diagonal edges (the
diagonal is always present and is not a model-selection choice).

**Criterion.** Each network has its own number of candidate parameters:
$p_\mu = d^2$, $p_\Delta = Kd^2$, $p_M = d(d-1)/2$, $p_E = Kd(d-1)/2$. Then

$$
\boxed{\;\text{EBIC}_\gamma = -2\hat\ell + \text{df}\,\log N + 2\gamma\Big[\text{df}_\mu\log p_\mu + \text{df}_\Delta\log p_\Delta + \text{df}_M\log p_M + \text{df}_E\log p_E\Big]\;}
$$

with $N = \sum_k N_k$ and $\gamma \in [0, 1]$ (`gamma_ebic`, default 0.5).
$\gamma = 0$ is the ordinary BIC; larger $\gamma$ penalizes model size more
heavily and is appropriate when the number of candidate edges is large relative
to $N$. The extra term is the log of the number of models of the given size
(up to constants), which is what makes EBIC consistent for support recovery in
the high-dimensional regime.

**Validity guard.** If any $S_{\varepsilon,k}(\hat B_k)$ is not positive
definite, the Gaussian likelihood is unbounded and any finite value of $\hat\ell$
is meaningless; such a model is assigned $\text{EBIC} = +\infty$ so it loses
every comparison, and its parameters are not used to warm-start the next grid
point. If every grid point is invalid the fit stops with an error advising
`pd_approx = TRUE`.

Two known second-order caveats: the likelihood is evaluated at the shrunken
estimates rather than at refitted maximum-likelihood values on the selected
support, which slightly favors smaller penalties; and an edge whose value is
split between $M$ and one or more $E_k$ is counted once per component.

---

## 12. Searching the four-dimensional penalty lattice

The selection is fully *coupled*: every evaluated candidate is the complete
joint model, obtained by running Algorithm A to convergence at that penalty
quadruple and scoring it with the joint EBIC. (An earlier staged design that
selected the temporal penalties without reference to the precision network was
rejected because it decouples the two halves of the model.)

Let the grids have sizes $n_\mu, n_\Delta, n_M, n_E$ (defaults 20, 20, 8, 8;
the two precision axes share `n_lambda_prec`). Two search strategies are
offered (`search`):

**`"grid"`** evaluates all $n_\mu n_\Delta n_M n_E$ lattice points in nested
loops (outermost $\lambda_\mu$, innermost $\lambda_E$) and returns the global
lattice minimizer. Each point is warm-started from the previously evaluated
point, temporal parameters and ADMM state alike (a *continuation* path), and
only from points with finite EBIC. The cost is the product of the four grid
sizes, so the defaults ($20 \cdot 20 \cdot 8 \cdot 8 = 25{,}600$ fits) are too
many in practice; a message is printed above 1000 points. With small
per-axis sizes (e.g. $5^4 = 625$ or $8^4 = 4096$) the search is feasible.

**`"coordinate"`** alternates two-dimensional sweeps:

> Start at $(\lambda_\mu, \lambda_\Delta)$ = top of both temporal grids and $(\lambda_M, \lambda_E)$ = midpoints of the precision grids.
> Repeat (at most `max_sweeps`, default 10):
> 1. Evaluate all $n_\mu n_\Delta$ temporal pairs at the current precision pair; move to the face minimizer.
> 2. Evaluate all $n_M n_E$ precision pairs at the new temporal pair; move to the face minimizer.
> 3. Stop when the selected quadruple no longer changes.

Every evaluation is still the full joint model; only the walk through the
lattice differs. Evaluated points are cached, so revisits are free. The result
is a blockwise (coordinate-wise) optimum on the lattice; since face minima are
monotone, the final point is the best of all visited points. The cost is a few
times $n_\mu n_\Delta + n_M n_E$, i.e. the *sum* rather than the product of the
face sizes. In an $8^4$ benchmark this reduced 4096 fits to 252 (a factor of
about 20) and selected the identical quadruple. Unvisited entries of the
returned `ebic_grid` are `NA`.

---

## 13. Special cases and relatives

**Ordinary multitask VAR (`gvar = FALSE`).** Setting $\Omega_k = I$ for all $k$
and dropping the precision block turns Algorithm A into a single call of the
temporal solver: the objective becomes the multitask least-squares problem

$$
\sum_k \frac{N_k}{2}\big\|\text{residuals}_k\big\|^2 + \lambda_\mu\|\mu\|_1 + \lambda_\Delta\sum_k\|\Delta_k\|_1
$$

(in covariance form, $\sum_k [\tfrac12 \operatorname{tr}(B_k\mathsf{G}_kB_k^\top) - \operatorname{tr}(B_k\mathsf{H}_k^\top)]$), which is the common-plus-individual VAR of Fisher et al. (2022) fitted by the same FISTA solver with the $\Omega$ factor switched off. Because there is no precision matrix in that objective, its EBIC uses a single *pooled* innovation covariance $\hat\Sigma = \frac{1}{N}\sum_k N_k S_{\varepsilon,k}(\hat B_k)$ and the profiled log-likelihood $\hat\ell = -\frac{N}{2}(\log\det\hat\Sigma + d)$, with only the two temporal df terms. The penalties are selected on the 2-D $(\lambda_\mu, \lambda_\Delta)$ grid with warm starts along $\lambda_\Delta$ and row seeding along $\lambda_\mu$. The `gvar = TRUE` machinery is a strict generalization: at $\Omega_k = I$ the weighted gradient, objective and Lipschitz operator reduce entrywise to the unweighted ones.

**Zero penalties.** At $\lambda_\mu = \lambda_\Delta = 0$ the temporal
subproblem is unpenalized and, whatever $\Omega_k$ is, the gradient
$N_k\Omega_k(B_kS_{0,k} - S_{10,k})$ vanishes at the per-subject Yule–Walker
solution. The decomposition into $\mu$ and $\Delta_k$ is then arbitrary
(Section 3.2), but the $B_k$ are the single-subject estimates. This is a useful
unit test and is what the package's tests check.

**One subject ($K = 1$).** With $K = 1$ and $\lambda_\Delta$ large enough to
zero $\Delta_1$, and $\lambda_E$ large enough to zero the off-diagonal of
$E_1$, the model collapses to a single-subject graphical VAR with transition
$\mu$ and precision $M + \operatorname{diag}(E_1)$; the alternation is then a
covariance-only analogue of the MRCE algorithm (Rothman, Levina and Zhu, 2010),
which is what `fit_graphical_var()` implements for one subject.

**Gaussian data.** With all-Gaussian, standardized data the link is the
identity, $S_0$ and $S_{10}$ are ordinary sample autocovariances, and the whole
procedure is a penalized Gaussian VAR with shared structure across subjects.
The copula layer adds nothing and removes nothing; it is only when some
marginals are discrete that Section 2 does real work.

---

## 14. The complete algorithm

Putting the pieces together, `fit_multitask(obj, gvar = TRUE)` does the
following.

> **Input.** $K$ data matrices ($n_k \times d$), marginal families, `pd_approx`, penalty type, grid sizes or supplied penalties, $\gamma$, `search`, $\rho$, tolerances.
>
> **Stage 0: per-subject latent moments** (constructor `timecop_multitask()`).
> For each subject: estimate marginal parameters; compute observed lag-0 and lag-1 cross-covariances; build the link coefficients from Hermite expansions; invert the link entrywise by spline interpolation; clamp to $[-1, 1]$; optionally repair the joint $2d \times 2d$ latent covariance to be PD. Store $S_{0,k}$, $S_{10,k}$, $N_k = n_k - 1$.
>
> **Stage 1: weights.** If `penalty = "adaptive"`, compute the per-subject pilot and the four weight matrices (Section 10); otherwise all weights are one.
>
> **Stage 2: grids.** Form $\mathsf{G}_k = N_kS_{0,k}$, $\mathsf{H}_k = N_kS_{10,k}$. Build the four log-spaced penalty sequences from their (weight-aware) KKT anchors (Section 9), using the Yule–Walker residual covariances for the precision anchors. Supplied penalties override.
>
> **Stage 3: sanity check.** Warn if any subject's joint latent covariance is not PSD.
>
> **Stage 4: search.** For each visited lattice point $(\lambda_\mu, \lambda_\Delta, \lambda_M, \lambda_E)$, in grid or coordinate order:
> * run Algorithm A (alternation), warm-started from the previous valid point:
>   * temporal block: FISTA with monotone restart on the stacked $(\mu, \Delta_1, \dots, \Delta_K)$, gradient $N_k\Omega_k(B_kS_{0,k} - S_{10,k})$, Lipschitz constant by power iteration;
>   * precision block: ADMM with the eigenvalue $\Omega$-update, edgewise soft-threshold $(M, E)$-update, and scaled dual update;
>   * until the temporal parameters stop moving;
> * score by the joint EBIC (Section 11);
> * keep the running minimizer.
>
> **Stage 5: report.** Return the EBIC-minimizing $\hat\mu, \hat\Delta_k, \hat B_k, \hat\Omega_k, \hat M, \hat E_k$, the four selected penalties and sequences, the EBIC surface, the four df counts, the number of outer iterations, and the final objective value; message if a selection sits on the boundary of an automatic grid.

---

## 15. Map from mathematics to code

| Mathematical object | Where it lives |
|---|---|
| $S_{0,k}$, $S_{10,k}$ | `subjects[[k]]@cov_z_hat[,, p + 1]` and `[,, p]` of the `timecop_multitask` object |
| $N_k$ | slot `N` (`n_k - p`) |
| Link $L_{ij}$, Hermite coefficients $g_{i,m}$ | `latent_var_link()`, `link_coefs()`, `hermite_coefs()`; `Polys[[m]]` stores $\mathrm{He}_{m-1}$ |
| Inverse link $L_{ij}^{-1}$ | `latent_var_invlink()` via `interpolation()` and `nat_spline()` |
| PD repair of $J$ | `check_pd()` (alternating projections) |
| $S_{\varepsilon,k}(B)$ | `gvar_resid_cov(B, S0, S0, t(S10), S10)` |
| $\mathsf{G}_k, \mathsf{H}_k$ | `G[[k]] = N[k] * S0[[k]]`, `H[[k]] = N[k] * t(S10[[k]])` in `fit_multitask()` |
| Temporal solver (Algorithm B) | `multitask_pgd()` |
| $\Omega$-update (Section 8.3) | `multitask_omega_update()` |
| $(M, E)$-update (Section 8.4) | `multitask_common_unique_update()` |
| Precision ADMM (Section 8) | `multitask_glasso_admm()` |
| Outer alternation (Algorithm A) | `multitask_gvar_alternate()` |
| KKT anchors and grids (Section 9) | `multitask_prec_anchors()`, `multitask_lambda_grids()` |
| Adaptive weights (Section 10) | `multitask_adaptive_weights()`, `gvar_weights_A()`, `gvar_weights_Omega()` |
| Joint EBIC (Section 11) | `multitask_gvar_ebic()` |
| Pooled EBIC for `gvar = FALSE` | `multitask_ebic()` |
| Grid / coordinate search (Section 12) | the `gvar` branch of `fit_multitask()` (`eval_point()` closure) |

**Orientation convention.** The temporal solver works with transposed
parameter blocks. Its stacked working variable holds $\mu^\top$ and
$\Delta_k^\top$, and correspondingly it consumes $\mathsf{H}_k^\top = N_kS_{10,k}^\top$
(the `H` list in the code). In this transposed space the gradient block is
$(\mathsf{G}_k\,\Theta_k^\top - \mathsf{H}_k^\top)\,\Omega_k$, i.e. the
precision matrix appears as a *right* multiplication; transposing recovers
$\Omega_k(\Theta_k\mathsf{G}_k - \mathsf{H}_k)$ of Section 7.2 because
$\mathsf{G}_k$ and $\Omega_k$ are symmetric. Every external interface of the
solver (warm starts, weights, outputs) is in natural orientation; the transposes
happen on entry and exit. The only place this leaks is in the grid
construction, where the temporal gradients are divided by $W^\top$ rather than
$W$.

**Objects returned by `fit_multitask(gvar = TRUE)`** (class
`timecop_multitask_gvar_fit`): `mu_hat`, `delta_hat` (list), `B_hat` (list),
`Omega_hat` (list, PD), `M_hat`, `E_hat` (list), `lambda_mu`, `lambda_delta`,
`lambda_M`, `lambda_E`, `ebic`, `ebic_grid` (4-D array), the four
`lambda_*_seq`, `df` (named vector), `gamma_ebic`, `penalty`, `search`,
`outer_iter`, `obj` (final value of $F$), and `object` (the input).

**Minimal usage.**

```r
library(timecop)

# K subjects, each an n_k x d matrix; one family per variable, shared across subjects
obj <- timecop_multitask(
  data      = list(X1, X2, X3),
  family    = list("Bernoulli", "Poisson", "Gaussian"),
  pd_approx = TRUE
)

fit <- fit_multitask(
  obj,
  gvar           = TRUE,
  penalty        = "adaptive",
  n_lambda_mu    = 8, n_lambda_delta = 8, n_lambda_prec = 8,
  search         = "coordinate",
  verbose        = TRUE
)

fit$mu_hat        # common temporal network
fit$delta_hat[[2]]  # subject 2's temporal deviation
fit$M_hat         # common contemporaneous edges (zero diagonal)
fit$E_hat[[2]]    # subject 2's unique edges + precision diagonal
fit$df            # nonzero counts in the four networks
```

---

## 16. References

* Beck, A. and Teboulle, M. (2009). A fast iterative shrinkage-thresholding algorithm for linear inverse problems. *SIAM Journal on Imaging Sciences*, 2(1), 183–202.
* Boyd, S., Parikh, N., Chu, E., Peleato, B. and Eckstein, J. (2011). Distributed optimization and statistical learning via the alternating direction method of multipliers. *Foundations and Trends in Machine Learning*, 3(1), 1–122.
* Chen, J. and Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. *Biometrika*, 95(3), 759–771.
* Danaher, P., Wang, P. and Witten, D. M. (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. *Journal of the Royal Statistical Society: Series B*, 76(2), 373–397.
* Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties. *Journal of the American Statistical Association*, 96(456), 1348–1360.
* Fisher, Z. F., Kim, Y., Fredrickson, B. L. and Pipiras, V. (2022). Penalized estimation and forecasting of multiple subject intensive longitudinal data. *Psychometrika*, 87, 1–29.
* Foygel, R. and Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. *Advances in Neural Information Processing Systems*, 23.
* Friedman, J., Hastie, T. and Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. *Biostatistics*, 9(3), 432–441.
* Jia, Y., Kechagias, S., Livsey, J., Lund, R. and Pipiras, V. (2023). Latent Gaussian count time series. *Journal of the American Statistical Association*, 118(541), 596–606.
* O'Donoghue, B. and Candès, E. (2015). Adaptive restart for accelerated gradient schemes. *Foundations of Computational Mathematics*, 15, 715–732.
* Rothman, A. J., Levina, E. and Zhu, J. (2010). Sparse multivariate regression with covariance estimation. *Journal of Computational and Graphical Statistics*, 19(4), 947–962.
* Zou, H. (2006). The adaptive lasso and its oracle properties. *Journal of the American Statistical Association*, 101(476), 1418–1429.
