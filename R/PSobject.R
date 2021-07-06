PSobject <- function(S.formula, Y.formula, Y.family, monotonicity = "default", ER = c(), trunc = FALSE){
  symbol_list <- c()
  param_list <- list()
  S_list <- list()
  Y_list <- list()
  strata <- c()
  .count <- 1
  
  if (monotonicity == "Default" || monotonicity == "default")
    strata <- c("00", "01", "11")
  else if (monotonicity == "None" || monotonicity == "none")
    strata <- c("00", "01", "10", "11")
  else if (monotonicity == "Strong" || monotonicity == "strong")
    strata <- c("00", "01")
  else
    stop("monotonicity option must be one of default, none and strong.")
  
  if (trunc){
    if ("00" %in% ER || "01" %in% ER || "10" %in% ER)
      stop("Exclusion restriction can only hold for stratum 11 under truncation.")
  }
  
  for (stratum in ER){
    if (!(stratum %in% strata))
      warning(paste0("Stratum ", stratum, " does not exist by monotonicity settings. ER is ignored for this stratum."))
  }
  
  .getRHS <- function(formula){
    return (formula[[length(formula)]])
  }
  
  .getLHS <- function(formula){
    stopifnot(length(formula) == 3)
    return (formula[[2]])
  }
  
  .extractRHS <- function(formula){
    update(formula, 1 ~ .)
  }
  
  .extractLHS <- function(formula){
    update(formula, . ~ 1)
  }
  
  LHS_terms <- function(formula){
    all.vars(attr(terms(.extractLHS(formula)), "variables"))
  }
  
  ## deal with outcome, treatment, intervention varnames.
  treatment_var <- LHS_terms(S.formula)[1]
  intervention_var <- LHS_terms(S.formula)[2]
  outcome_var <- LHS_terms(Y.formula)
  
  ## deal with S.formula
  
  S.RHS <- .extractRHS(S.formula)
  for (stratum in strata){
    .tmp <- expand_formula(S.RHS, symbol_list, .count)
    tmp <- .tmp$formula
    symbol_list <- .tmp$symbol_list
    S_list[[stratum]] <- list(
      name = as.character(tmp),
      fun = get_mean(as.character(tmp))
    )
    if (attr(tmp, "intercept")){
      param_list[[.count]] <- list(
        name = gen_name('S', stratum, 'Intercept'),
        model = "normal", 
        args = c("mean" = 0, "sd" = 100))
      .count <- .count + 1
    }
    if (attr(tmp, "count") != 0){
      for (i in .count:(.count + attr(tmp, "count") - 1)){
        param_list[[i]] <- list(
          name = gen_name('S', stratum, 'Coef', attr(tmp, "term.names")[i + 1 - .count]),
          model = "normal", 
          args = c("mean" = 0, "sd" = 10))
      }
      .count <- .count + attr(tmp, "count")
    }
  }
  
  Y.RHS <- .extractRHS(Y.formula)
  for (stratum in strata){
    for (z in 0:1){
      if (trunc){
        if (substring(stratum, z + 1, z + 1) == '0'){
          Y_list[[stratum]][[1 + z]] <- list(
            name = "0",
            fun = function(x, p) NA,
            invlink = function(x) x
          )
          next
        }
      }
      .tmp <- expand_formula(Y.RHS, symbol_list, .count)
      tmp <- .tmp$formula
      symbol_list <- .tmp$symbol_list
      eta <- as.character(tmp)
      if (attr(tmp, "intercept")){
        param_list[[.count]] <- list(
          name = gen_name('Y', c(stratum, z), 'Intercept'),
          model = "normal", 
          args = c("mean" = 0, "sd" = 100))
        .count <- .count + 1
      }
      if (attr(tmp, "count") != 0){
        for (i in .count:(.count + attr(tmp, "count") - 1)){
          param_list[[i]] <- list(
            name = gen_name('Y', c(stratum, z), 'Coef', attr(tmp, "term.names")[i + 1 - .count]),
            model = "normal", 
            args = c("mean" = 0, "sd" = 10))
        }
        .count <- .count + attr(tmp, "count")
      }
      if (Y.family$family == "gaussian") {
        sigma <- paste0('@', .count)
        param_list[[.count]] <- list(
          name = gen_name('Y', c(stratum, z), 'Sigma'),
          model = "inv_gamma", args = c("alpha" = 0.1, "beta" = 0.1))
        .count <- .count + 1
        if (Y.family$link == "identity") {
          str <- eta
        }
        else if (Y.family$link == "log"){
          str <- paste0("exp(", eta, ")")
        }
        else if (Y.family$link == "inverse"){
          str <- paste0("1 / (", eta, ")")
        }
        else
          stop("Link function must be one of identity, log and inverse.")
        Y_list[[stratum]][[1 + z]] <- list(
          name = paste0("normal_lpdf(. | ", str, ", ", sigma, ")"),
          fun = get_mean(eta),
          invlink = Y.family$linkinv
        )
      }
      else if (Y.family$family == "binomial"){
        if (Y.family$link == "logit") {
          str <- paste0("inv_logit(", eta, ")")
        }
        else if (Y.family$link == "probit"){
          str <- paste0("inv_Phi(", eta, ")")
        }
        else if (Y.family$link == "cauchit"){
          str <- paste0("atan(", eta, ") / pi() + .5")
        }
        else if (Y.family$link == "log"){
          str <- paste0("exp(", eta, ")")
        }
        else if (Y.family$link == "cloglog"){
          str <- paste0("inv_cloglog(", eta, ")")
        }
        else
          stop("Link function must be one of logit, probit, cauchit, log and cloglog.")
        Y_list[[stratum]][[1 + z]] <- list(
          name = paste0("bernoulli_lpmf(. | ", str, ")"),
          fun = get_mean(eta),
          invlink = Y.family$linkinv
        )
      }
      else if (Y.family$family == "Gamma"){
        shape <- paste0('@', .count)
        param_list[[.count]] <- list(
          name = gen_name('Y', c(stratum, z), 'Alpha'),
          model = "inv_gamma", 
          args = c("alpha" = 0.1, "beta" = 0.1))
        .count <- .count + 1
        if (Y.family$link == "identity") {
          str <- paste0(shape, " / (", eta, ")")
        }
        else if (Y.family$link == "log"){
          str <- paste0(shape, " / exp(", eta, ")")
        }
        else if (Y.family$link == "inverse"){
          str <- paste0(shape, " * (", eta, ")")
        }
        else
          stop("Link function must be one of identity, log and inverse.")
        Y_list[[stratum]][[1 + z]] <- list(
          name = paste0("gamma_lpdf(. | ", shape, ", ", str, ")"),
          fun = get_mean(eta),
          invlink = Y.family$linkinv
        )
      }
      else if (Y.family$family == "poisson"){
        if (Y.family$link == "identity") {
          str <- eta
        }
        else if (Y.family$link == "log"){
          str <- paste0("exp(", eta, ")")
        }
        else if (Y.family$link == "sqrt"){
          str <- paste0("(", eta, ")^2")
        }
        else
          stop("Link function must be one of identity, log and sqrt.")
        Y_list[[stratum]][[1 + z]] <- list(
          name = paste0("poisson_lpmf(. | ", str, ")"),
          fun = get_mean(eta),
          invlink = Y.family$linkinv
        )
      }
      else if (Y.family$family == "inverse.gaussian"){
        lambda <- paste0('@', .count)
        param_list[[.count]] <- list(
          name = gen_name('Y', c(stratum, z), 'Lambda'),
          model = "inv_gamma", 
          args = c("alpha" = 0.1, "beta" = 0.1))
        .count <- .count + 1
        if (Y.family$link == "identity") {
          str <- eta
        }
        else if (Y.family$link == "log"){
          str <- paste0("exp(", eta, ")")
        }
        else if (Y.family$link == "inverse"){
          str <- paste0("1 / (", eta, ")")
        }
        else if (Y.family$link == "1/mu^2"){
          str <- paste0("1 / (", eta, ")^2")
        }
        else
          stop("Link function must be one of identity, log, inverse and 1/mu^2.")
        Y_list[[stratum]][[1 + z]] <- list(
          name = paste0("inv_gaussian_lpdf(. | ", str, ", ", lambda, ")"),
          fun = get_mean(eta),
          invlink = Y.family$linkinv
        )
      }
      else
        stop("'Y.family' must be one of gaussian, binomial, Gamma, poisson and inverse.gaussian.")
      if (stratum %in% ER){
        Y_list[[stratum]][[2]] <- Y_list[[stratum]][[1]]
        break
      }
    }
  }
  
  if (Y.family$family %in% c("gaussian"))
    Y_type = "continuous"
  else if (Y.family$family %in% c("binomial"))
    Y_type = "binary"
  else if (Y.family$family %in% c("Gamma", "inverse.gaussian"))
    Y_type = "positive"
  else if (Y.family$family %in% c("poission"))
    Y_type = "count"
  
  result <- list(
    outcome_var = outcome_var,
    treatment_var = treatment_var,
    intervention_var = intervention_var,
    outcome_type = Y_type,
    symbol_list = symbol_list,
    param_list = param_list,
    S_list = S_list,
    Y_list = Y_list,
    strata = strata,
    ER = ER
  )
  class(result) <- "PSobject"
  return (result)
}

print.PSobject <- function(obj){
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

