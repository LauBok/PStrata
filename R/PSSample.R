PSSample <- function(filename, data, ...) {
  rstan::stan(file = filename, data = data, ...)
}