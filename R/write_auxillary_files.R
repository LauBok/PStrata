write.auxillary.files <- function() {
  ### create folder 
  if (!dir.exists("auto_generated_files"))
    dir.create("auto_generated_files")
  
  ### family info
  fileConn <- file("auto_generated_files/family_info.txt")
  lines <- c(
    "gaussian            real        normal_lpdf         1   sigma   positive",
    "binomial            binary      bernoulli_lpmf      0",
    "Gamma               positive    Gamma_lpdf          1   alpha   positive",
    "poisson             count       poisson_lpmf        0",
    "inverse.gaussian    real        inv_gaussian_lpdf   1   lambda  real"
  )
  writeLines(lines, fileConn, sep = "\n")
  close(fileConn)
  
  ### function_table
  fileConn <- file("auto_generated_files/function_table.txt")
  lines <- c(
    "gaussian            identity    normal_lpdf         identity_func       1   sigma   positive",
    "gaussian            log         normal_lpdf         exp                 1   sigma   positive",
    "gaussian            inverse     normal_lpdf         inv_func            1   sigma   positive",
    "binomial            logit       bernoulli_lpmf      inv_logit           0",
    "binomial            probit      bernoulli_lpmf      inv_Phi             0",
    "binomial            cauchit     bernoulli_lpmf      inv_cauchit_func    0",
    "binomial            log         bernoulli_lpmf      exp                 0",
    "binomial            cloglog     bernoulli_lpmf      inv_cloglog         0",
    "Gamma               inverse     Gamma_lpdf          inv_func            1   alpha   positive",
    "Gamma               identity    Gamma_lpdf          identity_func       1   alpha   positive",
    "Gamma               log         Gamma_lpdf          exp                 1   alpha   positive",
    "poisson             log         poisson_lpmf        exp                 0",
    "poisson             identity    poisson_lpmf        identity_func       0",
    "poisson             sqrt        poisson_lpmf        square_func         0",
    "inverse.gaussian    1/mu^2      inv_gaussian_lpdf   inv_square_func     1   lambda  real",
    "inverse.gaussian    inverse     inv_gaussian_lpdf   inv_func            1   lambda  real",
    "inverse.gaussian    identity    inv_gaussian_lpdf   identity_func       1   lambda  real",
    "inverse.gaussian    log         inv_gaussian_lpdf   exp                 1   lambda  real"
  )
  writeLines(lines, fileConn, sep = "\n")
  
  ### link_info
  fileConn <- file("auto_generated_files/link_info.txt")
  lines <- c(
    "gaussian            identity    identity_func",
    "gaussian            log         exp",
    "gaussian            inverse     inv_func",
    "binomial            logit       inv_logit",
    "binomial            probit      inv_Phi",
    "binomial            cauchit     inv_cauchit_func",
    "binomial            log         exp",
    "binomial            cloglog     inv_cloglog",
    "Gamma               inverse     inv_func",
    "Gamma               identity    identity_func",
    "Gamma               log         exp",
    "poisson             log         exp",
    "poisson             identity    identity_func",
    "poisson             sqrt        square_func",
    "inverse.gaussian    1/mu^2      inv_square_func",
    "inverse.gaussian    inverse     inv_func",
    "inverse.gaussian    identity    identity_func",
    "inverse.gaussian    log         exp"
  )
  writeLines(lines, fileConn, sep = "\n")
  
  ### link_info
  fileConn <- file("auto_generated_files/function_implement.txt")
  lines <- c(
    "<<< identity_func >>>",
    "real identity_func(real x) {",
    "    return x;",
    "}",
    "<<<>>>",
    "",
    "<<< inv_func >>>",
    "real inv_func(real x) {",
    "    return 1 / x;",
    "}",
    "<<<>>>",
    "",
    "<<< square_func >>>",
    "real square_func(real x) {",
    "    return x^2;",
    "}",
    "<<<>>>",
    "",
    "<<< inv_square_func >>>",
    "real inv_square_func(real x) {",
    "    return 1 / x^2;",
    "}",
    "<<<>>>",
    "",
    "<<< inv_cauchit_func >>>",
    "real inv_cauchit_func(real x) {",
    "    return atan(x) / pi() + 0.5;",
    "}",
    "<<<>>>",
    "",
    "<<< Gamma_lpdf >>>",
    "real Gamma_lpdf(real x, real mu, real alpha) {",
    "    return gamma_lpdf(x | alpha, alpha / mu);",
    "}",
    "<<<>>>",
    "",
    "<<< inv_gaussian_lpdf >>>",
    "real inv_gaussian_lpdf(real x, real mu, real lambda) {",
    "    real constant = log(lambda) / 2.0 - log(2 * pi()) / 2.0;",
    "    real kernel = -1.5 * log(x) - lambda * pow(x - mu, 2) / (2 * x * pow(mu, 2));",
    "    return constant + kernel;",
    "}",
    "<<<>>>"
  )
  writeLines(lines, fileConn, sep = "\n")
}