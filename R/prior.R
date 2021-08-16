prior_uniform <- function() {
  return (
    list(
      name = "uniform",
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