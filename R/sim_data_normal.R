#' Simulated Dataset for Normal Outcome
#'
#' A dataset generated for illustration of the principal stratification analysis.
#' This dataset represents the common case of non-compliance.
#' 
#' @details 
#' The dataset represents the scenario where actual treatment might not be in compliance
#' with the randomized (assigned) treatment. Defiers are ruled out, leaving three
#' strata, "never taker", "complier" and "always taker" randomly sampled with 
#' probability 0.3, 0.2 and 0.5 respectively. The assigned treatment \eqn{Z} is randomized
#' with 0.5 probability for either arm. The outcome \eqn{Y} is given by the following.
#' \describe{
#' \item{never taker}{\eqn{Y \sim N(3, 1)}}
#' \item{complier}{\eqn{Y \sim N(-1-Z, 0.5)}}
#' \item{always taker}{\eqn{Y \sim N(1, 2)}}
#' }
#' The exclusion restriction assumption holds for never takers and always takers in this
#' generated dataset.
#'
#' @format ## `sim_data_normal`
#' A data frame with 1,000 rows and 4 columns:
#' \describe{
#'   \item{S}{Principal Strata: "never taker", "complier" or "always taker"}
#'   \item{Z}{Randomized treatment arm: 0 = control, 1 = treatment}
#'   \item{D}{Actual treatment arm: 0 = control, 1 = treatment}
#'   \item{Y}{Outcome}
#' }
"sim_data_normal"