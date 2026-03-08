#' Prior functions
#'
#' Define prior functions used in \code{PStrata}.
#'
#' @param mu,sigma,df,alpha,beta parameters for the prior distribution
#'
#' @return A list, including the following items.
#' \describe{
#' \item{name}{name of the distribution}
#' \item{type}{type of the distribution, one character string of "real" or "positive"}
#' \item{args}{a named list of all the input parameters}
#' \item{call}{a function call object of the prior distribution on the parameters}
#' }
#'
#' @name prior
#' @rdname prior
NULL

#' Construct a prior distribution object (internal factory)
#' @param name Stan distribution name
#' @param type "real" or "positive"
#' @param ... Named arguments in Stan argument order
#' @return A prior list
#' @noRd
make_prior <- function(name, type, ...) {
  args <- list(...)
  call_obj <- if (name == "flat") {
    NULL
  } else {
    as.call(c(list(as.name(name)), unname(args)))
  }
  list(name = name, type = type, args = args, call = call_obj)
}

#' @rdname prior
#' @export
prior_flat <- function() {
  make_prior("flat", "real")
}

#' @rdname prior
#' @export
prior_normal <- function(mu = 0, sigma = 1) {
  make_prior("normal", "real", mu = mu, sigma = sigma)
}

#' @rdname prior
#' @export
prior_t <- function(mu = 0, sigma = 1, df = 1) {
  make_prior("t", "real", df = df, mu = mu, sigma = sigma)
}

#' @rdname prior
#' @export
prior_cauchy <- function(mu = 0, sigma = 1) {
  make_prior("cauchy", "real", mu = mu, sigma = sigma)
}

#' @rdname prior
#' @export
prior_lasso <- function(mu = 0, sigma = 1) {
  make_prior("double_exponential", "real", mu = mu, sigma = sigma)
}

#' @rdname prior
#' @export
prior_logistic <- function(mu = 0, sigma = 1) {
  make_prior("logistic", "real", mu = mu, sigma = sigma)
}

#' @rdname prior
#' @export
prior_chisq <- function(df = 1) {
  make_prior("chi_square", "positive", df = df)
}

#' @rdname prior
#' @export
prior_inv_chisq <- function(df = 1) {
  make_prior("inv_chi_square", "positive", df = df)
}

#' @rdname prior
#' @export
prior_exponential <- function(beta = 1) {
  make_prior("exponential", "positive", beta = beta)
}

#' @rdname prior
#' @export
prior_gamma <- function(alpha = 1, beta = 1) {
  make_prior("gamma", "positive", alpha = alpha, beta = beta)
}

#' @rdname prior
#' @export
prior_inv_gamma <- function(alpha = 1, beta = 1) {
  make_prior("inv_gamma", "positive", alpha = alpha, beta = beta)
}

#' @rdname prior
#' @export
prior_weibull <- function(alpha = 1, sigma = 1) {
  make_prior("weibull", "positive", alpha = alpha, sigma = sigma)
}
