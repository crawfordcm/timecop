## Script to generate example dataset for timecop package
library(timecop)

set.seed(42)
sim <- latent_var_sim(
  d = 2, n = 200, p = 1,
  param = list(0.5, 0.5),
  phi_lv = matrix(c(0.4, 0.2, 0.2, 0.4), 2, 2),
  family = list("Bernoulli", "Bernoulli")
)

timecop_example <- t(sim$X_t)  # 200 x 2 matrix

usethis::use_data(timecop_example, overwrite = TRUE)
