PSObject <- function(
  S.model, Y.model, Y.family, monotonicity, ER,
  prior_intercept = prior_uniform(),
  prior_coefficient = prior_normal(),
  prior_sigma = prior_inv_gamma(),
  prior_alpha = prior_inv_gamma(),
  prior_lambda = prior_inv_gamma(),
  prior_theta = prior_normal()
) {
  
  get_list_vars <- function(S.model, Y.model){
    symbol_S <- manipulate_formula(S.model)$symbols$RHS
    symbol_Y <- manipulate_formula(Y.model)$symbols$RHS
    return (unique(c(symbol_S, symbol_Y)))
  }
  
  # covariates
  symbol_list <- get_list_vars(S.model, Y.model)
  
  # analyze monotonicity
  if (monotonicity == "default")
    strata <- c("00", "01", "11")
  else if (monotonicity == "none")
    strata <- c("00", "01", "10", "11")
  else if (monotonicity == "strong")
    strata <- c("00", "01")
  
  # create model for each stratum
  parameter_list <- list()
  stratum_model_list <- list()
  outcome_model_list <- list()
  start_num <- 1
  
  # models for stratum
  S_model_info <- manipulate_formula(S.model)
  treatment.var <- S_model_info$symbols$LHS[[1]]
  intervention.var <- S_model_info$symbols$LHS[[2]]
  
  for (stratum in strata) {
    if (stratum == strata[1]){
      stratum_model_list[[stratum]] <- 0
    }
    else{
      new_param <- S_model_info$parameters
      for (i in 1:length(new_param)){
        new_param[[i]]$name <- paste(
          'S', stratum, new_param[[i]]$name, sep = '_'
        )
      }
      parameter_list <- c(parameter_list, new_param)
      stratum_model_list[[stratum]] <- S_model_info$evaluate(start_num)
      start_num <- start_num + S_model_info$num_of_parameters
    }
  }
  
  #models for outcome
  Y_model_info <- manipulate_formula(Y.model)
  outcome.var <- Y_model_info$symbols$LHS[[1]]
  if (F) # to change
    censor.var <- Y_model_info$symbols$LHS[[2]]
  else
    censor.var <- NULL
  
  for (stratum in strata) {
    model <- model_template(
      Y.model, Y.family, stratum %in% ER, stratum
    )(start_num)
    parameter_list <- c(parameter_list, model$param_list)
    outcome_model_list[[stratum]] <- model$model_list
    start_num <- model$start
  }
  
  # outcome type
  if (Y.family$family %in% c("gaussian"))
    Y_type = "continuous"
  else if (Y.family$family %in% c("binomial"))
    Y_type = "binary"
  else if (Y.family$family %in% c("Gamma", "inverse.gaussian"))
    Y_type = "positive"
  else if (Y.family$family %in% c("poission"))
    Y_type = "count"
  else if (Y.family$family %in% c("survival"))
    Y_type = "survival"
  
  for (i in 1:length(parameter_list)) {
    prior_type <- parameter_list[[i]]$prior_type
    prior_dist <- eval(str2lang(prior_type))
    parameter_list[[i]]$prior_dist <- prior_dist
  }
  
  return (structure(list(
    variables = list(
      treatment = treatment.var,
      intervention = intervention.var,
      outcome = outcome.var,
      censor = censor.var
    ),
    symbol_list = symbol_list,
    parameter_list = parameter_list,
    stratum_model_list = stratum_model_list,
    outcome_model_list = outcome_model_list,
    outcome_type = Y_type,
    strata = strata
  ), class = "PSObject"))
}

write.pso <- function(obj, filename = NULL){
  substitute_par <- function(call_object, env = c()){
    return (eval(substitute(
      substitute(y, env), 
      list(y = call_object))
    ))
  }
  
  unname_symbol <- function(call_object, env_vars) {
    return (substitute_par(call_object, env_vars))
  }
  
  reformat <- function(call_object){
    args <- sapply(call_object[-1], deparse, backtick = F)
    str_args <- paste0(args, collapse = ', ')
    strings <- c(
      as.character(call_object[[1]]),
      '(. | ',
      str_args,
      ')'
    )
    return (paste0(strings, collapse = ''))
  }
  
  if (!is.null(filename))
    fileConn <- file(filename)
  
  env <- lapply(
    1:length(obj$symbol_list), 
    function(i) str2lang(paste('`$', i, '`', sep = ''))
  )
  names(env) <- obj$symbol_list
  
  indent <- function(str, indent) {
    return (paste0(paste0(rep(' ', indent), collapse = ''), str))
  }
  lines <- c()
  # Y type
  lines <- c(lines, paste0("Y: ", obj$outcome_type))
  # S info
  lines <- c(lines, paste0("S: ", paste0(strtoi(obj$strata, base = 2), collapse = ' ')))
  # Covariates
  lines <- c(lines, "covariate {")
  if (length(obj$symbol_list) > 0){
    for (i in 1:length(obj$symbol_list)){
      lines <- c(
        lines, 
        indent(paste0('$', i, ': ', obj$symbol_list[i]), indent = 4)
      )
    }
  }
  lines <- c(lines, "}")
  # Parameters
  lines <- c(lines, "parameter {")
  for (i in 1:length(obj$parameter_list)){
    lines <- c(
      lines, 
      indent(
        paste0('@', i, ': <', 
               obj$parameter_list[[i]]$type, '> ',
               obj$parameter_list[[i]]$name
        ), indent = 4
      )
    )
  }
  lines <- c(lines, "}")
  # Prior Distribution
  lines <- c(lines, "prior {")
  for (i in 1:length(obj$parameter_list)){
    call <- obj$parameter_list[[i]]$prior_dist$call
    if (is.null(call))
      tmp_str <- "uniform()"
    else
      tmp_str <- deparse(call)
    lines <- c(
      lines, 
      indent(paste0('@', i, ' ~ ', tmp_str), indent = 4)
    )
  }
  lines <- c(lines, "}")
  # Strata models
  lines <- c(lines, "strata {")
  for (i in 1:length(obj$stratum_model_list)){
    lines <- c(
      lines, 
      indent(
        paste0(
          strtoi(names(obj$stratum_model_list)[i], 2), 
          ": ", 
          deparse(
            unname_symbol(obj$stratum_model_list[[i]], env),
            backtick = F
          )
        ), indent = 4
      )
    )
  }
  lines <- c(lines, "}")
  # Outcome models
  lines <- c(lines, "outcome {")
  for (i in 1:length(obj$outcome_model_list)){
    for (j in 0:1){
      lines <- c(
        lines, 
        indent(
          paste0(
            strtoi(names(obj$outcome_model)[i], 2), 
            ", ", j,  ": ", 
            reformat(
              unname_symbol(
                obj$outcome_model[[i]][[j + 1]]$model, env
              )
            )
          ), indent = 4
        )
      )
    }
  }
  lines <- c(lines, "}")
  
  if (!is.null(filename)) {
    writeLines(lines, fileConn)
    close(fileConn)
  }
  return (lines)
}

print_future.PSObject <- function(obj){
  width = 80
  cat(nice.print("Principle Stratification Object", width, 3, '*'))
  cat(nice.print(paste0("Treatment : ", obj$treatment_var, '   ', 
                        "Intervention : ", obj$intervention_var, '   ', 
                        "Outcome : ", obj$outcome_var, 
                        ' (', obj$outcome_type,')'), width, 2, ' '))
  cat('\n')
  cat(nice.print("Covariate Names", width, 3, '-'))
  cur_col_num = 5
  cat("     ")
  if (length(obj$symbol_list) != 0){
    for (i in 1:length(obj$symbol_list)){
      this.length = stringr::str_length(i) + 3 + stringr::str_length(obj$symbol_list[i])
      if (cur_col_num + this.length >= 85){
        cat("\n     ")
        cur_col_num = 5
      }
      cat(paste0('$', i, ": ", obj$symbol_list[i]))
      cur_col_num = cur_col_num + this.length
      if (cur_col_num < 25){
        cat(paste0(rep(' ', 25 - cur_col_num), collapse = ''))
        cur_col_num = 25
      }
      else if (cur_col_num < 45){
        cat(paste0(rep(' ', 45 - cur_col_num), collapse = ''))
        cur_col_num = 45
      }
      else if (cur_col_num < 65){
        cat(paste0(rep(' ', 65 - cur_col_num), collapse = ''))
        cur_col_num = 65
      }
      else{
        cur_col_num = 85
      }
    }
  }  
  cat('\n')
  cat(nice.print('', width, 0))
  cat('\n')
  
  cat(nice.print("Parameter Names", width, 3, '-'))
  cur_col_num = 5
  cat("     ")
  for (i in 1:length(obj$param_list)){
    this.length = stringr::str_length(i) + 3 + stringr::str_length(obj$param_list[[i]]$name)
    if (cur_col_num + this.length >= 85){
      cat("\n     ")
      cur_col_num = 5
    }
    cat(paste0('@', i, ": ", obj$param_list[[i]]$name))
    cur_col_num = cur_col_num + this.length
    if (cur_col_num < 25){
      cat(paste0(rep(' ', 25 - cur_col_num), collapse = ''))
      cur_col_num = 25
    }
    else if (cur_col_num < 45){
      cat(paste0(rep(' ', 45 - cur_col_num), collapse = ''))
      cur_col_num = 45
    }
    else if (cur_col_num < 65){
      cat(paste0(rep(' ', 65 - cur_col_num), collapse = ''))
      cur_col_num = 65
    }
    else{
      cur_col_num = 85
    }
  }
  cat('\n')
  cat(nice.print('', width, 0))
  cat('\n')
  
  cat(nice.print("Parameter List", width, 3, '-'))
  cur_col_num = 5
  cat("     ")
  for (i in 1:length(obj$param_list)){
    this.length = stringr::str_length(i) + 3 + stringr::str_length(obj$param_list[[i]]$model) + 2 + 
      sum(stringr::str_length(obj$param_list[[i]]$args)) + length(obj$param_list[[i]]$args) + 
      sum(stringr::str_length(attr(obj$param_list[[i]]$args, "names"))) + 2 * (length(obj$param_list[[i]]$args) - 1)
    if (cur_col_num + this.length >= 85){
      cat("\n     ")
      cur_col_num = 5
    }
    cat(paste0('@', i, ": ", obj$param_list[[i]]$model, '('))
    cat(paste0(paste(
      attr(obj$param_list[[i]]$args, "names"),
      obj$param_list[[i]]$args,
      sep = '='
    ), collapse = ", "))
    cat(')')
    cur_col_num = cur_col_num + this.length
    if (cur_col_num < 45){
      cat(paste0(rep(' ', 45 - cur_col_num), collapse = ''))
      cur_col_num = 45
    }
    else{
      cur_col_num = 85
    }
  }
  cat('\n')
  cat(nice.print('', width, 0))
  cat('\n')
  cat(nice.print('Strata', width, 3))
  n_keep <- stringr::str_length(obj$treatment_var) + 8
  half_width <- (width - 6 - n_keep) %/% 2
  cat(paste0(rep(' ', n_keep), collapse = ''))
  cat('++')
  cat(nice.print(paste0(obj$treatment_var, "[1] = 0"), half_width, 2, '+', newline = F))
  cat('++')
  cat(nice.print(paste0(obj$treatment_var, "[1] = 1"), half_width, 2, '+', newline = F))
  cat('++')
  cat('\n')
  cat(paste0(paste0(rep(' ', n_keep), collapse = ''), '++', 
             paste0(rep(' ', half_width), collapse = ''), '++', 
             paste0(rep(' ', half_width), collapse = ''), '++', '\n', collapse = ''))
  cat(paste0(paste0(obj$treatment_var, "[0] = 0 "), '++', 
             nice.print(paste0(
               ifelse("00" %in% obj$strata, "Never-taker", ""),
               ifelse('00' %in% obj$ER, '*', '')), half_width, 0, ' ', newline = F), '++', 
             nice.print(paste0(
               ifelse("01" %in% obj$strata, "Complier", ""),
               ifelse('01' %in% obj$ER, '*', '')), half_width, 0, ' ', newline = F), '++',
             '\n', collapse = ''))
  cat(paste0(paste0(rep(' ', n_keep), collapse = ''), '++', 
             paste0(rep(' ', half_width), collapse = ''), '++', 
             paste0(rep(' ', half_width), collapse = ''), '++', '\n', collapse = ''))
  cat(paste0(paste0(rep(' ', n_keep), collapse = ''), 
             paste0(rep('+', half_width * 2 + 6), collapse = ''), '\n'))
  cat(paste0(paste0(rep(' ', n_keep), collapse = ''), '++', 
             paste0(rep(' ', half_width), collapse = ''), '++', 
             paste0(rep(' ', half_width), collapse = ''), '++', '\n', collapse = ''))
  cat(paste0(paste0(obj$treatment_var, "[0] = 1 "), '++', 
             nice.print(paste0(
               ifelse("10" %in% obj$strata, "Defier", ""),
               ifelse('10' %in% obj$ER, '*', '')), half_width, 0, ' ', newline = F), '++', 
             nice.print(paste0(
               ifelse("11" %in% obj$strata, "Always-taker", ""),
               ifelse('11' %in% obj$ER, '*', '')), half_width, 0, ' ', newline = F), '++',
             '\n', collapse = ''))
  cat(paste0(paste0(rep(' ', n_keep), collapse = ''), '++', 
             paste0(rep(' ', half_width), collapse = ''), '++', 
             paste0(rep(' ', half_width), collapse = ''), '++', '\n', collapse = ''))
  cat(paste0(paste0(rep(' ', n_keep), collapse = ''), 
             paste0(rep('+', half_width * 2 + 6), collapse = ''), '\n'))
  cat(paste0(paste0(rep(' ', n_keep), collapse = ''), ' * Exclusion Restriction holds.\n'))
  
  cat('\n')
  cat(nice.print("S-model", width))
  for (i in 1:length(obj$S_list)){
    cat(paste0('log Pr(S = ', names(obj$S_list)[i], ') = ', obj$S_list[[i]]$name, ' + C\n'))
  }
  cat('\n')
  
  cat(nice.print("Y-model", width))
  for (i in 1:length(obj$S_list)){
    cat(paste0('log Pr(. | S = ', names(obj$Y_list)[i], ', Z = 0) = ', obj$Y_list[[i]][[1]]$name, ' + C\n'))
    cat(paste0('log Pr(. | S = ', names(obj$Y_list)[i], ', Z = 1) = ', obj$Y_list[[i]][[2]]$name, ' + C\n'))
  }
  cat('\n')
}