#' Extract important information from an R formula object
#' 
#' Parse the formula to obtain the names of the responses and covariates, and the design matrix
#' with \code{data} provided. When the formula includes random effects, specified using the syntax
#' for package \code{lme4}, the design matrix corresponding for each cluster variable is also provided.
#' @param formula a formula, syntaxed as \code{lm()} or \code{lme4::lmer()}.
#' @param data a data frame, including variables specified in \code{formula}.
#' @return A list, containing at least the following components:
#' \describe{
#' \item{formula}{formula as the input argument}
#' \item{non_random_formula}{formula with only fixed effects}
#' \item{random_effect}{boolean value of whether a random effect is included in the formula}
#' \item{response}{symbols that appear on the left hand side of the formula}
#' \item{has_intercept}{boolean value of whether intercept is included in the formula}
#' \item{predictors}{symbols of all the linear predictors}
#' \item{num_of_predictors}{number of linear predictors}
#' \item{model_matrix}{the design matrix given by the fixed effects of the model}
#' }
#' @return When random effects are included, the following item is also included:
#' \describe{
#' \item{random_effect_list}{a list of all the random effects, each of which is a list with a design matrix for the linear terms specific to the effect, the name of these linear terms, and the levels that the random effect can be}
#' }
#' @seealso \code{something}
parse.formula <- function(formula, data) {
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
  terms <- terms.formula(non_random_formula)
  has_intercept <- attr(terms, "intercept")
  model_matrix <- model.matrix(non_random_formula, data)
  
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
    formula = formula, 
    non_random_formula = non_random_formula,
    response = LHS_symbol,
    random_effect = random_effect,
    has_intercept = has_intercept,
    predictors = colnames(model_matrix),
    num_of_predictors = ncol(model_matrix),
    model_matrix = model_matrix,
    random_matrix_list = random_matrix_list
  ))
}
