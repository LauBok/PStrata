parse.formula <- function(formula, data) {
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
  model_matrix <- model.matrix(
    update.formula(formula, ~ . + 1), 
    dplyr::mutate_if(data, function(x) is.numeric(x) && var(x) != 0, scale)
  )[, -1, drop = F]
  return (list(
    formula = formula, 
    response = LHS_symbol,
    has_intercept = has_intercept,
    predictors = colnames(model_matrix),
    num_of_predictors = ncol(model_matrix),
    model_matrix = model_matrix
  ))
}
