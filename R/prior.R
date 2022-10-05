#' Prior functions
#' 
#' Define prior functions used in \code{PStrata}.
#' 
#' @export prior_flat
#' @export prior_normal
#' @export prior_t
#' @export prior_cauchy
#' @export prior_lasso
#' @export prior_logistic
#' @export prior_chisq
#' @export prior_inv_chisq
#' @export prior_exponential
#' @export prior_gamma
#' @export prior_inv_gamma
#' @export prior_weibull
#' @param ... parameters for the prior distribution
#' @return A list, including the following items.
#' \describe{
#' \item{name}{name of the distribution}
#' \item{type}{type of the distribution, one character string of "real" or "positive"}
#' \item{args}{a named list of all the input parameters}
#' \item{call}{a function call object of the prior distribution on the parameters}
#' }
#' 

prior_flat <- function() {
  return (
    list(
      name = "flat",
      type = "real",
      args = list(),
      call = NULL
    )
  )
}

prior_normal <- function(mu = 0, sigma = 1){
  return (
    list(
      name = "normal",
      type = "real", 
      args = list(mu = mu, sigma = sigma),
      call = as.call(list(quote(normal), mu, sigma))
    )
  )
}

prior_t <- function(mu = 0, sigma = 1, df = 1){
  return (
    list(
      name = "t",
      type = "real", 
      args = list(mu = mu, sigma = sigma, df = df),
      call = as.call(list(quote(t), df, mu, sigma))
    )
  )
}

prior_cauchy <- function(mu = 0, sigma = 1){
  return (
    list(
      name = "cauchy",
      type = "real", 
      args = list(mu = mu, sigma = sigma),
      call = as.call(list(quote(cauchy), mu, sigma))
    )
  )
}

prior_lasso <- function(mu = 0, sigma = 1) {
  return (
    list(
      name = "double_exponential",
      type = "real", 
      args = list(mu = mu, sigma = sigma),
      call = as.call(list(quote(double_exponential), mu, sigma))
    )
  )
}

prior_logistic <- function(mu = 0, sigma = 1) {
  return (
    list(
      name = "logistic",
      type = "real", 
      args = list(mu = mu, sigma = sigma),
      call = as.call(list(quote(logistic), mu, sigma))
    )
  )
}

prior_chisq <- function(df = 1) {
  return (
    list(
      name = "chi_square",
      type = "positive", 
      args = list(df = df),
      call = as.call(list(quote(chi_square), df))
    )
  )
}

prior_inv_chisq <- function(df = 1) {
  return (
    list(
      name = "inv_chi_square",
      type = "positive", 
      args = list(df = df),
      call = as.call(list(quote(inv_chi_square), df))
    )
  )
}

prior_exponential <- function(beta = 1) {
  return (
    list(
      name = "exponential",
      type = "positive", 
      args = list(beta = beta),
      call = as.call(list(quote(exponential), beta))
    )
  )
}

prior_gamma <- function(alpha = 1, beta = 1) {
  return (
    list(
      name = "gamma",
      type = "positive", 
      args = list(alpha = alpha, beta = beta),
      call = as.call(list(quote(gamma), alpha, beta))
    )
  )
}

prior_inv_gamma <- function(alpha = 1, beta = 1) {
  return (
    list(
      name = "inv_gamma",
      type = "positive", 
      args = list(alpha = alpha, beta = beta),
      call = as.call(list(quote(inv_gamma), alpha, beta))
    )
  )
}

prior_weibull <- function(alpha = 1, sigma = 1) {
  return (
    list(
      name = "weibull",
      type = "positive", 
      args = list(alpha = alpha, sigma = sigma),
      call = as.call(list(quote(weibull), alpha, sigma))
    )
  )
}
