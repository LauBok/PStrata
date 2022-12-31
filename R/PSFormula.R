#' Set up a model formula for use in \pkg{PStrata}
#' 
#' Set up a model formula for use in \pkg{PStrata} package allowing users to specify
#' the treatment indicator, the post-randomization confounding variables, the outcome variable, and possibly the covariates.
#' For survival outcome, a censoring indicator is also specified.
#' Users can also define (potentially non-linear) transforms of the covariates and include random effects for clusters.
#' 
#' 
#' @param formula an object of class \code{formula} (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted. The details of model specification are given in 'Details'.
#' @param data a data frame containing the variables named in \code{formula}.
#' 
#' @details 
#' Two models are required for the principal stratification analysis: the principal stratum model and the outcome model.
#' 
#' \subsection{General formula structure}{
#' For the principal stratum model, the \code{formula} argument accepts formulas of the following syntax:
#' 
#' \code{treatment + postrand ~ terms}
#' 
#' The \code{treatment} variable refers to the name of the binary treatment indicator.
#' The \code{postrand} variable refers to the name of the binary post-randomization confounding variable.
#' The \code{terms} part includes all of the predictors used for the principal stratum model.
#' 
#' For the outcome model, the \code{formula} argument accepts formulas of the similar syntax:
#' 
#' \code{response [+ observed] ~ terms}
#' 
#' The \code{response} variable refers to the name of the outcome variable.
#' The \code{terms} part includes all of the predictors used for the outcome model. 
#' The \code{observed} variable shall not be used for ordinary response.
#' When the true response is subject to right censoring (also called survival outcome in relevant literature),
#' the \code{response} variable should refer to the observed or censored response, and the \code{observed} variable should
#' be an indicator of whether the true response is observed. 
#' For example, suppose the true time for an event is \eqn{T} and the time of censoring is \eqn{C},
#' Then, the \code{response} variable should refer to \eqn{\min(T, C)}, the actual time of the event or censoring, whichever comes earlier,
#' and the indicator \code{observed} is 1 if \eqn{T < C} and 0 otherwise.
#' 
#' The \code{terms} specified in the principal stratum model and the outcome model can be different.
#' }
#' 
#' \subsection{Multiple post-randomization confounding variables}{
#' If multiple post-randomization confounding variables exist, one can specify all of them using 
#' the following syntax:
#' 
#' \code{treatment + postrand_1 + postrand_2 + ... + postrand_n ~ terms}
#' 
#' The post-randomization confounding variables are provided in place of \code{postrand_1} to
#' \code{postrand_n}. Up to this version, all of these variables should be binary indicators.
#' Note that the order of these post-randomization confounding variables will not
#' affect the result of the estimation of the parameters, but it will be important
#' in specifying other parameters, such as \code{strata} and \code{ER} (see \code{\link{PStrata}}).
#' }
#' 
#' \subsection{Non-linear transformation of the predictors}{
#' The syntax for the predictors follow the conventions as used in \code{link{formula}}.
#' The part \code{terms} consists of a series of terms concatenated by \code{+}, 
#' each term being the name of a variable, or the interaction of several variables separated by \code{:}.
#' 
#' Apart from \code{+} and \code{:}, a number of other operators are also useful.
#' The \code{*} operator is a short-hand for factor crossing: 
#' \code{a*b} is interpreted as \code{a + b + a:b}.
#' The \code{^} operator means factor crossing to a specific degree. For example,
#' \code{(a + b + c)^2} is interpreted as \code{(a + b + c) * (a + b + c)},
#' which is identical to \code{a + b + c + a:b + a:c + b:c}.
#' The \code{-} operator removes specified terms, so that \code{(a + b + c)^2 - a:b} is
#' identical to \code{a + b + c + a:c + b:c}.
#' The \code{-} operator can be also used to remove the intercept term, such as 
#' \code{x - 1}. One can also use \code{x + 0} to remove the intercept term.
#' 
#' Arithmetic expressions such as \code{a + log(b)} are also legal.
#' However, arithmetic expressions may contain special symbols that are defined for other use, such as \code{+}, \code{*}, \code{^} and \code{-}.
#' To avoid confusion, the function \code{\link{I}()} can be used to bracket portions where the operators should be interpreted in arithmetic sense.
#' For example, in \code{x + I(y + z)}, the term \code{y + z} is interpreted as the sum of \code{y} and \code{z}.
#' }
#' 
#' \subsection{Group level random effect}{
#' When effects assumed to vary across grouping variables are considered, one can 
#' specify such effects by adding terms in the form of \code{gterms | group}, where
#' \code{group} refers to the group indicator (usually a \code{factor}), and
#' \code{gterms} specifies the terms whose coefficients are group-specific, drawn 
#' from a population normal distribution. 
#' 
#' The most common situation for group level random effect is to include group-specific
#' intercepts to account for unmeasured confounding. 
#' For example, \code{x + y + (1 | g)} specifies a model with population predictors
#' \code{x} and \code{y}, as well as random intercept for each level of \code{g}.
#' 
#' For more complex random effect structures, refer to \code{\link[lme4:lmer]{lme4::lmer}}.
#' However, structures other than simple random intercepts and slopes may lead to unexpected behaviors.
#' }
#' 
#' @returns \code{PSFormula} returns an object of class \code{PSFormula}, 
#' which is a \code{list} containing for following components.
#' \describe{
#'  \item{\code{full_formula}}{input formula as is}
#'  \item{\code{data}}{input data frame}
#'  \item{\code{fixed_eff_formula}}{input formula with only fixed effects}
#'  \item{\code{response_names}}{character vector with names of variables that appear on the left hand side of input formula}
#'  \item{\code{has_random_effect}}{logical indicating whether random effects are specified in the input formula}
#'  \item{\code{has_intercept}}{logical indicating whether the input formula has an intercept}
#'  \item{\code{fixed_eff_names}}{character vector with names of all variables included as fixed effects}
#'  \item{\code{fixed_eff_count}}{integer indicating the number of variables (factors are converted to and counted as dummy variables)}
#'  \item{\code{fixed_eff_matrix}}{fixed-effect design matrix}
#'  \item{\code{random_eff_list}}{a list containing information for each random effect. 
#'  Such information is a list with the corresponding design matrix, the term names and the factor levels.}
#' }
#' 
#' @examples 
#' df <- data.frame(
#'   X = 1:10, 
#'   Z = c(0,0,0,0,0,1,1,1,1,1),
#'   D = c(0,0,0,1,1,1,0,0,1,1),
#'   R = c(1,1,1,1,2,2,2,3,3,3)
#'  )
#' PSFormula(Z + D ~ X + I(X^2) + (1 | R), df)
#' 
#' @seealso 
#' \code{\link{formula}}, \code{\link[lme4:lmer]{lmer}}.
#' 
#' @export
PSFormula <- function(formula, data) {
  symbols_AST <- function(AST) {
    if (is.name(AST))
      return (as.character(AST))
    else if (is.call(AST))
      return (unlist(unique(sapply(AST[-1], symbols_AST))))
    else 
      return (NULL)
  }
  
  if (is.null(lme4::findbars(formula))) {
    random_effect <- F
  }
  else {
    random_effect <- T
  }
  
  non_random_formula <- lme4::nobars(formula)
  
  LHS <- if(length(formula) == 3) formula[[2]] else NULL
  LHS_symbol <- symbols_AST(LHS)
  terms <- stats::terms.formula(non_random_formula)
  has_intercept <- as.logical(attr(terms, "intercept"))
  model_matrix <- stats::model.matrix(non_random_formula, data)
  
  random_matrix_list <- NULL
  if (random_effect) {
    lf <- lme4::lFormula(formula, data)
    random_matrix_list <- mapply(
      function(x, y, z) list(
        matrix = t(as.matrix(x)),
        terms = y,
        factors = z
      ),
      lf$reTrms$Ztlist, 
      lf$reTrms$cnms,
      lf$reTrms$flist,
      SIMPLIFY = F
    )
  }
  
  return (list(
    full_formula = formula, 
    data = data,
    fixed_eff_formula = non_random_formula,
    response_names = LHS_symbol,
    has_random_effect = random_effect,
    has_intercept = has_intercept,
    fixed_eff_names = colnames(model_matrix),
    fixed_eff_count = ncol(model_matrix),
    fixed_eff_matrix = model_matrix,
    random_eff_list = random_matrix_list
  ))
}