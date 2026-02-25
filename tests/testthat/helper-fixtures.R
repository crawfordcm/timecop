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
