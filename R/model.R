survival <- function(link = "identity") {
  return (list(
    family = "survival",
    link = link
  ))
}

model <- function(family, formula) {
  family_chr <- family$family
  link_chr <- family$link
  formula <- formula
  return (list(
    family = family_chr,
    link = link_chr,
    formula = formula
  ))
}

get_param_from_model <- function(model) {
  prse_fml <- parse.formula(model$formula)
  param_list <- list()
  if (prse_fml$has_intercept) {
    param_list <- c(param_list, list(list(
      type = "real", prior_type = "prior_intercept", name = "Intrcpt"
    )))
  }
  param_list <- c(param_list, list(list(
    type = "real_vct", dim = prse_fml$num_of_predictors,
    prior_type = "prior_coefficient", name = "Coef"
  )))
  if (model$family == "gaussian")
    param_list <- c(param_list, list(list(
      type = "positive", prior_type = "prior_sigma", name = "Sigma"
    )))
  else if (model$family == "Gamma")
    param_list <- c(param_list, list(list(
      type = "positive", prior_type = "prior_alpha", name = "Alpha"
    )))
  else if (model$family == "inv.gaussian")
    param_list <- c(param_list, list(list(
      type = "positive", prior_type = "prior_lambda", name = "Lambda"
    )))
  else if (model$family == "survival")
    param_list <- c(param_list, list(list(
      type = "real", prior_type = "prior_theta", name = "Theta"
    )))
  return (param_list)
}

get_dist_str <- function(family, Y.formula, group){
  family_name <- family$family
  link <- family$link
  prse_fml_Y <- parse.formula(Y.formula)
  name_intrcpt <- paste(c('Y', group, "Intrcpt"), collapse = '_')
  name_coef <- paste(c('Y', group, "Coef"), collapse = '_')
  str_core <- c()
  if (prse_fml_Y$has_intercept)
    str_core <- c(str_core, name_intrcpt)
  str_core <- c(str_core, paste0("$ * ", name_coef))
  str_mean <- paste(str_core, collapse = ' + ')
  group_prefix <- paste(c("Y", group, ""), collapse = '_')
  
  if (family_name == "gaussian"){
    str_sigma <- paste0(group_prefix, "Sigma")
    if (link == "identity")
      str_inner <- str_mean
    else if (link == 'qlog')
      str_inner <- paste0("exp(", str_mean, ")")
    else if (link == 'inverse')
      str_inner <- paste0("1/(", str_mean, ")")
    return (paste0("normal_lpdf(. | ", str_inner, ", ", str_sigma, ")"))
  }
  if (family_name == "binomial") {
    if (link == 'logit')
      str_inner <- paste0("inv_logit(", str_mean, ")")
    else if (link == 'probit')
      str_inner <- paste0("inv_Phi(", str_mean, ")")
    else if (link == "cauchit")
      str_inner <- paste0("atan(", str_mean, ") / pi + 0.5")
    else if (link == "log")
      str_inner <- paste0("exp(", str_mean, ")")
    else if (link == "cloglog")
      str_inner <- paste0("inv_cloglog(", str_mean, ")")
    return (paste0("bernoulli_lpmf(. | ", str_inner, ")"))
  }
  if (family_name == "Gamma") {
    str_alpha <- paste0(group_prefix, "Alpha")
    if (link == "identity")
      str_inner <- paste0(str_alpha, " / (", str_mean, ")")
    else if (link == "log")
      str_inner <- paste0(str_alpha, "/exp(", str_mean, ")")
    else if (link == "inverse")
      str_inner <- paste0(str_alpha, " * (", str_mean, ")")
    return (paste0("gamma_lpdf(. | ", str_alpha, ", ", str_inner, ")"))
  }
  if (family_name == "poisson") {
    if (link == "identity")
      str_inner <- paste0(str_mean)
    else if (link == "log")
      str_inner <- paste0("exp(", str_mean, ")")
    else if (link == "sqrt")
      str_inner <- paste0("(", str_mean, ") ^ 2")
    return (paste0("poisson_lpmf(. | ", str_inner, ")"))
  }
  if (family_name == "inv.gaussian") {
    str_lambda <- paste0(group_prefix, "Lambda")
    if (link == "identity")
      str_inner <- paste0(str_mean)
    else if (link == "log")
      str_inner <- paste0("exp(", str_mean, ")")
    else if (link == "inverse")
      str_inner <- paste0("1 / (", str_mean, ")")
    else if (link == "1/mu^2")
      str_inner <- paste0("1 / (", str_mean, ")^2")
    return (paste0("inv_gaussian_lpdf(. | ", str_inner, ", ", str_lambda, ")"))
  }
  if (family_name == "survival") {
    str_theta <- paste0(group_prefix, "Theta")
    str_inner <- paste0(str_theta, " + ", str_mean)
    return (paste0("survival_lpdf(. | ", str_inner, ", ", str_theta, ", .)"))
  }
}
