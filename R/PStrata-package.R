#' PStrata: Principal STRATification Analysis for Data with Post-Randomization Confounding
#'
#' The \bold{PStrata} package is designed for estimating causal effects in the presence of post-treatment confounding using principal stratification.
#' It provides an interface to fit the Bayesian principal stratification model, which is a complex mixture model, using \bold{Stan}, a C++ package for obtaining full Bayesian inference.
#' The formula syntax is an extended version of the syntax applied in many regression functions and packages, such as \code{\link{lm}}, \code{\link{glm}} and \pkg{\link[lme4]{lme4-package}}, to provide a simple interface.
#' A wide variety of distributions and link functions are supported, allowing users to fit linear, binary or count data, and survival models with principal stratification.
#' Further modeling options include multiple post-treatment confounding variables and cluster random effects.
#' The monotonicity and exclusion restriction assumptions can be easily applied, and
#' prior specifications are flexible and encourage users to reflect their prior belief.
#' In addition, all parameters can be inferred from the posterior distribution, which enables further analysis other than provided by the package.
#'
#' @name PStrata_PACKAGE
#' @details
#' The Bayesian principal stratification analysis relies on two models, the principal stratum model and the outcome model.
#' The main workflow is:
#' \enumerate{
#'   \item \code{\link{PStrataModel}}: specify the model (formulas, strata, priors).
#'   \item \code{\link{fit}}: compile and run MCMC sampling via Stan.
#'   \item \code{\link{estimate}}: extract posterior potential outcomes.
#'   \item \code{\link{contrast}}: compute causal effect contrasts.
#' }
#' Based on the supplied formulas, data and additional information,
#' it automatically generates the Stan code via \code{\link{make_stancode}} and fits the model using \pkg{Stan}.
#'
#' The estimated probability for each principal stratum and the estimated mean response are calculated with \pkg{Stan} as it is faster and more space-efficient.
#' Post-processing methods \code{\link{summary}} and \code{\link{plot}} provide overviews and visualization.
#'
#' Because \pkg{PStrata} heavily relies on \pkg{Stan} for posterior sampling, a C++ compiler is required.
#' The program \pkg{Rtools} (available on \url{https://cran.r-project.org/bin/windows/Rtools/}) comes with a C++ compiler for Windows.
#' On Mac, Xcode is suggested. For further instructions on how to get the compilers running, please refer to the prerequisites section at the \href{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}{RStan-Getting-Started} page.
#'
#' @references
#' The Stan Development Team. Stan Modeling Language User's Guide and Reference Manual. \url{https://mc-stan.org/users/documentation/}
#'
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. \url{https://mc-stan.org/}
#' @seealso \code{\link{PStrataModel}}, \code{\link{fit}}, \code{\link{estimate}}, \code{\link{contrast}}
#' @importFrom stats coef formula nobs
NULL

# Suppress R CMD check NOTEs for ggplot2 aes() variables
utils::globalVariables(c("value", "Z", "S", "T", "mean", "label",
                         "2.5%", "97.5%"))
