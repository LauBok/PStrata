#' PStrata: Principal STRATification Analysis for Data with Post-Randomization Confounding
#' 
#' The \bold{PStrata} package is designed for estimating causal effects in the presense of post-treatment confounding using principal stratification.
#' It provides an interface to fit the Bayesian principal stratification model, which is a complex mixture model, using \bold{Stan}, a C++ package for obtaining full Bayesian inference.
#' The formula syntax is an extended version of the syntax applied in many regression functions and packages, such as \code{\link{lm}}, \code{\link{glm}} and \pkg{\link{lme4}}, to provide a simple interface. 
#' A wide variety of distributions and link functions are supported, allowing users to fit linear, binary or count data, and survival models with principal stratification.
#' Further modeling options include multiple post-treatment confounding variables and cluster random effects.
#' The monotonicity and exclusion restriction assumptions can be easily applied, and
#' prior specifications are flexible and encourage users to reflect their prior belief.
#' In addition, all parameters can be inferred from the posterior distribution, which enables further analysis other than provided by the package.
#' A frequentist weighting-based triply-robust estimator is also implemented for both ordinary outcomes and survival outcomes.
#' 
#' @name PStrata-package
#' @docType package
#' @alias PStrata-package
#' @useDynLib PStrata
#' @details 
#' The Bayesian principal stratification analysis relies on two models, the principal stratum model and the outcome model.
#' The main function of \pkg{PStrata} is \code{\link{PStrata}}, which uses formula syntax to specify these models.
#' Based on the supplied formulas, data and additional information allowing users to specify assumptions and prior distributions,
#' it automatically generates the Stan code via \code{\link{make_stancode}} and \code{\link{make_standata}}, and fits the model using \pkg{Stan}.
#' 
#' The estimated probability for each principal stratum and the estimated mean response are calculated with \pkg{Stan} as it is faster and more space-efficient. 
#' However, a large number of post-processing methods can also be applied.
#' \code{\link[=PStrata]{summary}} is perfectly suited for an overview of the estimated parameters, 
#' and \code{\link[=PStrata]{plot}} provides visualization of the principal stratification and the outcome distribution.
#' Customized outcome calculation can be performed with \link{TODO}.
#' 
#' Because \pkg{PStrata} heavily relies on \pkg{Stan} for posterior sampling, a C++ compiler is required.
#' The program \pkg{Rtools} (available on \url{https://cran.r-project.org/bin/windows/Rtools/}) comes with a C++ compiler for Windows.
#' On Mac, Xcode is suggested. For further instructions on how to get the compilers running, please refer to the prerequisites section at the \href{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}{RStan-Getting-Started} page.
#' 
#' @references 
#' The Stan Development Team. Stan Modeling Language User's Guide and Reference Manual. \url{https://mc-stan.org/users/documentation/}
#' 
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. \url{https://mc-stan.org/}
#' @seealso \code{\link{PStrata}}
NULL