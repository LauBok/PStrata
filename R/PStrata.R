#' Principal strata analysis for data with post-randomization intervention
#' 
#' description
#' 
#' @param S.formula a formula for the principal strata model, in the form of \code{Z + D1 + D2 ~ X1 + X2}, where \code{Z} is the treatment assignment and \code{D1} and \code{D2} are the post-randomization variables.
#' @param Y.formula a formula for the outcome model, in the form of \code{Y ~ X1 + X2}, where \code{Y} is the outcome variable
#' @param Y.family a family object, indicating the family of the model as used in \code{lm()}.
#' @param data a data frame object, the data used to do inference on.
#' @param strata a vector of integers, each integer referring to a stratum that is included in the model
#' @param ER a vector of integers, indicating the strata on which the exclusion restriction assumption is assumed.
#' @param prior_... the prior distribution for the corresponding parameters.
#' @return A text file, readable by \code{to_stan()}.

PStrata <- function(
    S.formula, Y.formula, Y.family, data, strata, ER,
    prior_intercept = prior_flat(),
    prior_coefficient = prior_normal(),
    prior_sigma = prior_inv_gamma(),
    prior_alpha = prior_inv_gamma(),
    prior_lambda = prior_inv_gamma(),
    prior_theta = prior_normal(),
    filename = "unnamed",
    ...
) {
  PSobject <- PSObject(S.formula, Y.formula, Y.family, data, strata, ER,
                       prior_intercept,
                       prior_coefficient,
                       prior_sigma,
                       prior_alpha,
                       prior_lambda,
                       prior_theta,
                       filename)
  to_stan(filename)
  post_samples <- PSSample(filename = paste0(filename, ".stan"),
           data = PSobject$stan_data,
           ...)
  res <- list(post_samples = post_samples)
  class(res) <- "PStrata"
  return (res)
}