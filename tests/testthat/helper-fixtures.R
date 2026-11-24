make_test_sim <- function(n = 200) {
  set.seed(42)
  latent_var_sim(
    d = 2, n = n, p = 1,
    param = list(0.5, 0.5),
    phi_lv = matrix(c(0.4, 0.2, 0.2, 0.4), 2, 2),
    family = list("Bernoulli", "Bernoulli")
  )
}

make_test_sim_poisson <- function(n = 200) {
  set.seed(42)
  latent_var_sim(
    d = 2, n = n, p = 1,
    param = list(3, 5),
    phi_lv = matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2),
    family = list("Poisson", "Poisson")
  )
}

# Small multi-subject (multitask) dataset: K Gaussian subjects sharing a common
# transition matrix, with an off-diagonal shared innovation precision edge.
make_multitask_data <- function(K = 2, d = 3, n = 300) {
  set.seed(7)
  mu <- matrix(0, d, d); diag(mu) <- 0.4
  if (d >= 2) mu[2, 1] <- 0.25
  M_true <- matrix(0, d, d)
  if (d >= 3) { M_true[1, 3] <- M_true[3, 1] <- -0.4 }
  Sig <- solve(diag(d) + M_true)
  fam <- as.list(rep("Gaussian", d))
  dat <- lapply(seq_len(K), function(k) {
    s <- latent_var_sim(d = d, n = n, p = 1, param = as.list(rep(0, d)),
                        phi_lv = mu, family = fam, sigma_lv = Sig)
    t(s$X_t)
  })
  list(data = dat, family = fam, d = d, K = K)
}
