#' Principal Stratification Analysis for Data with Post-Randomization Intervention
#' 
#' Perform pincipal stratification analysis when there are confounding variables
#' after randomization
#' 
#' @export PStrata
#' 
#' @param PSobject an object of class \code{\link{PSObject}}. 
#' If left blank, the object is constructed using the following arguments.
#' See \code{\link{PSObject}} for details.
#' @inheritParams PSObject
#' @param filename (optional) string. If not \code{NULL}, the stan file will be saved via
#' \code{\link{cat}} in a text file named after the string supplied.
#' @param ... additional parameters to be passed into \code{\link{PSSample}}.
#' 
#' @examples 
#' require(abind)
#' PSobj <- PSObject(
#'   S.formula = Z + D ~ 1,
#'   Y.formula = Y ~ 1,
#'   Y.family = gaussian("identity"),
#'   data = sim_data_normal,
#'   strata = c(n = "00*", c = "01", a = "11*")
#' )
#' 
#' #----- Not Run ------
#' # PStrata(PSobj, cores = 6, chains = 6, iter = 1000)
#' 
#' # Another example for survival data
#' PSobj <- PSObject(
#'   S.formula = Z + D ~ 1,
#'   Y.formula = Y + delta ~ 1,
#'   Y.family = survival("Cox"),
#'   data = sim_data_Cox,
#'   strata = c(`never-taker` = "00*", complier = "01")
#' )
#' 
#' #----- Not Run ------
#' # PStrata(PSobj, cores = 6, chains = 6, iter = 1000)
#' 
#' @return An object of class \code{PStrata} or \code{PStrata_survival}, 
#' which is a list containing the \code{PSObject} and the posterior samples (of class \code{stanfit})
PStrata <- function(
    PSobject,
    S.formula,
    Y.formula,
    Y.family,
    data = NULL,
    strata = NULL,
    ER = NULL,
    prior_intercept = prior_flat(),
    prior_coefficient = prior_normal(),
    prior_sigma = prior_inv_gamma(),
    prior_alpha = prior_inv_gamma(),
    prior_lambda = prior_inv_gamma(),
    prior_theta = prior_normal(),
    survival.time.points = 50,
    filename = NULL,
    ...
) {
  if (is.null(PSobject)){
    PSobject <- PSObject(S.formula, Y.formula, Y.family,
                         data, strata, ER, 
                         prior_intercept, 
                         prior_coefficient,
                         prior_alpha,
                         prior_lambda,
                         prior_theta,
                         survival.time.points)
  }
  stancode <- make_stancode(PSobject, filename)
  standata <- make_standata(PSobject)
  post_samples <- PSSample(
           data = standata,
           model_code = stancode,
           ...)
  res <- list(
    PSobject = PSobject,
    post_samples = post_samples
  )
  if (PSobject$Y.family$family == "survival_Cox" || PSobject$Y.family$family == "survival_AFT")
    class(res) <- "PStrata_survival"
  else
    class(res) <- "PStrata"
  return (res)
}

#' @export
print.PStrata <- function(x,...) {
  cat("PStrata Object with", x$PSobject$strata_info$num_strata, "strata.\n")
  cat("The estimated proportion for each strata: \n")
  strata_prob <- colMeans(rstan::extract(x$post_samples)$'strata_prob')
  names(strata_prob) <- x$PSobject$strata_info$strata_names
  cat(strata_prob, '\n')
  cat("Use summary() to show confidence intervals.\n")
  cat("Use PSOutcome() to show outcomes and PSContrast() to show constrasts.\n")
}

#' @export
print.PStrata_survival <- function(x, ...) {
  print.PStrata(x, ...)
}

#' @export
summary.PStrata <- function(object, ...) {
  strata_prob <- t(apply(rstan::extract(object$post_samples)$'strata_prob', 2,
                       function(x) c(mean(x), stats::sd(x), stats::quantile(x, 0.025), 
                                     stats::quantile(x, 0.25), stats::median(x), stats::quantile(x, 0.75),
                                     stats::quantile(x, 0.975))))
  rownames(strata_prob) <- object$PSobject$strata_info$strata_names
  colnames(strata_prob) <- c("mean", "sd", "2.5%", "25%", "median", "75%", "97.5%")
  return (strata_prob)
}

#' @export
summary.PStrata_survival <- function(object, ...) {
  return (summary.PStrata(object, ...))
}
