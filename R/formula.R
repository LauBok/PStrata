manipulate_formula <- function(formula){
  # Separate left and right side
  meta_formula <- function(formula) {
    len <- length(formula) # 2: Only RHS, 3: LHS and RHS
    LHS <- if(len == 3) formula[[2]] else NULL
    RHS <- formula[[len]]
    return (list(LHS = LHS, RHS = RHS))
  }
  
  # Obtain all symbols on both sides
  symbols_AST <- function(AST) {
    if (is.name(AST))
      return (as.list(AST))
    else if (is.call(AST))
      return (unlist(unique(sapply(AST[-1], symbols_AST))))
    else 
      return (NULL)
  }
  
  # Obtain all terms on one side
  terms_formula <- function(formula) {
    modify_AST <- function(AST){
      if (is.call(AST)){
        if (AST[[1]] == quote(`:`)) AST[[1]] <- quote(`*`)
        else if (AST[[1]] == quote(`^`)) AST[[1]] <- quote(pow)
        else if (AST[[1]] == quote(I)) AST <- modify_AST(AST[[2]])
        
        AST[-1] <- lapply(AST[-1], modify_AST)
      }
      return (AST)
    }
    terms <- terms.formula(formula)
    labels <- attr(terms, "term.labels") # all terms as a character vector
    # add 1 manually if intercept exists
    if (attr(terms, "intercept") == 1)
      labels <- c(1, labels)
    # then we need to convert these terms to call objects 
    return (sapply(sapply(labels, str2lang), modify_AST))
  }
  
  
  # parameter list
  parameter_list <- function(formula){
    params <- list()
    if (attr(terms.formula(formula), "intercept")){
      params <- c(
        params, 
        list(list(type = "real", prior_type = "prior_intercept", 
                  name = "Intercept"))
      )
    }
    labels <- attr(terms.formula(formula), "term.labels")
    for (label in labels){
      params <- c(
        params, list(list(type = "real", prior_type = "prior_coefficient",
                          name = paste0("Coef_", tidy_name(label))))
      )
    }
    return (params)
  }
  
  tidy_name <- function(name = NA){
    .res <- name
    .res <- gsub("I\\(|\\)", "", .res)
    .res <- gsub('\\^', "pow", .res)
    .res <- gsub('\\*', "mult", .res)
    .res <- gsub('\\+', "add", .res)
    .res <- gsub('\\-', "minus", .res)
    return (gsub("[^0-9a-zA-z_]|\\^", '_', .res))
  }
  
  # evaluate formula
  evaluate_formula <- function(start_num = 1L, env = list()) {
    terms <- terms_formula(formula)
    numbers <- start_num : (start_num - 1 + length(terms))
    coefficients <- sapply(paste('`@', numbers, '`', sep = ''), str2lang)
    multiply <- function(x, y) as.call(list(quote(`*`), x, y))
    terms_with_coef <- unname(mapply(multiply, coefficients, terms))
    
    # accumulate these terms together using `+`.
    accumulate <- function(x) {
      if (length(x) == 1) return (x[[1]])
      if (length(x) == 2) return (as.call(list(quote(`+`), x[[1]], x[[2]])))
      return (as.call(list(
        quote(`+`),
        accumulate(x[-length(x)]),
        x[[length(x)]]
      )))
    }
    
    return (eval(substitute(
      substitute(y, env), 
      list(y = accumulate(terms_with_coef)))
    ))
  }
  
  # return values
  meta_result = meta_formula(formula)
  result <- list(
    meta = meta_result,
    symbols = list(
      LHS = symbols_AST(meta_result$LHS),
      RHS = symbols_AST(meta_result$RHS)
    ),
    terms = terms_formula(formula),
    num_of_parameters = length(terms_formula(formula)),
    intercept = attr(terms.formula(formula), "intercept"),
    parameters = parameter_list(formula),
    evaluate = evaluate_formula
  )
  return (result)
}