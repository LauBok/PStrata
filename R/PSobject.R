PSObject <- function(
  S.formula, Y.formula, Y.family, data, strata, ER,
  prior_intercept = prior_flat(),
  prior_coefficient = prior_normal(),
  prior_sigma = prior_inv_gamma(),
  prior_alpha = prior_inv_gamma(),
  prior_lambda = prior_inv_gamma(),
  prior_theta = prior_normal(),
  filename = NULL
) {
  pso <- write.txt(S.formula, Y.formula, Y.family, strata, ER,
            prior_intercept, prior_coefficient,
            prior_sigma, prior_alpha, prior_lambda, prior_theta,
            filename)
  stan_data <- get.stan.data(S.formula, Y.formula, data)
  return (list(pso_file = pso, stan_data = stan_data))
}
