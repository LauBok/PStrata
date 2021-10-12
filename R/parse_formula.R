parse.formula <- function(formula) {
  symbols_AST <- function(AST) {
    if (is.name(AST))
      return (as.character(AST))
    else if (is.call(AST))
      return (unlist(unique(sapply(AST[-1], symbols_AST))))
    else 
      return (NULL)
  }
  
  LHS <- if(length(formula) == 3) formula[[2]] else NULL
  LHS_symbol <- symbols_AST(LHS)
  terms <- terms.formula(formula)
  has_intercept <- attr(terms, "intercept")
  predictors <- attr(terms, "term.labels")
  num_of_predictors <- length(predictors)
  return (list(
    formula = formula, 
    response = LHS_symbol,
    has_intercept = has_intercept,
    predictors = predictors,
    num_of_predictors = num_of_predictors
  ))
}