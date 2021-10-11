PSObject <- function(
  S.model, Y.model, Y.family, monotonicity, ER,
  prior_intercept = prior_uniform(),
  prior_coefficient = prior_normal(),
  prior_sigma = prior_inv_gamma(),
  prior_alpha = prior_inv_gamma(),
  prior_lambda = prior_inv_gamma(),
  prior_theta = prior_normal()
) {
  # analyze monotonicity
  if (monotonicity == "default")
    strata <- c("00", "01", "11")
  else if (monotonicity == "none")
    strata <- c("00", "01", "10", "11")
  else if (monotonicity == "strong")
    strata <- c("00", "01")
  
  # create model for each stratum
  parameter_list <- list()
  
  # models for stratum
  prse_fml_S <- parse.formula(S.model)
  treatment.var <- prse_fml_S$response[1]
  intervention.var <- prse_fml_S$response[2]
  S_has_intercept <- prse_fml_S$has_intercept
  S_num_of_param <- prse_fml_S$num_of_predictors
  
  for (stratum in strata[-1]) {
    if (S_has_intercept)
      parameter_list <- append(parameter_list, list(list(
        type = "real", prior_type = "prior_intercept", name = "Intrcpt",
        group = stratum, class = 'S'
      )))
    parameter_list <- append(parameter_list, list(list(
      type = "real_vct", dim = S_num_of_param,
      prior_type = "prior_coefficient", name = "Coef",
      group = stratum, class = 'S'
    )))
  }
  
  # models for outcome
  prse_fml_Y <- parse.formula(Y.model)
  outcome.var <- prse_fml_Y$response[1]
  censor.var <- prse_fml_Y$response[2]
  Y_params <- get_param_from_model(model(Y.family, Y.model))
  
  for (stratum in strata) {
    for (z in 0:1) {
      if (stratum %in% ER && z == 1)
        break
      Y_tmp_params <- Y_params
      for (i in 1:length(Y_tmp_params)) {
        Y_tmp_params[[i]]$group <- c(stratum, z)
        Y_tmp_params[[i]]$class <- 'Y'
      }
      parameter_list <- append(parameter_list, Y_tmp_params)
    }
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
  
  
  return (structure(list(
    variables = list(
      treatment = treatment.var,
      intervention = intervention.var,
      outcome = outcome.var,
      censor = censor.var
    ),
    parameter_list = parameter_list,
    outcome_type = Y_type,
    strata = strata,
    ER = ER,
    S.model = S.model, 
    Y.model = Y.model, 
    Y.family = Y.family,
    priors = list(
      prior_intercept = prior_intercept,
      prior_coefficient = prior_coefficient,
      prior_sigma = prior_sigma,
      prior_alpha = prior_alpha,
      prior_lambda = prior_lambda,
      prior_theta = prior_theta
    ),
    S.dim = prse_fml_S$num_of_predictors,
    Y.dim = prse_fml_Y$num_of_predictors
  ), class = "PSObject"))
}

obj <- PSObject(
  S.model = Z + D ~ X1 * X2,
  Y.model = Y + C ~ X3 + X4,
  Y.family = survival(),
  monotonicity = "default",
  ER = c('00', '11')
)

write.pso <- function(obj, filename = NULL){
  if (!is.null(filename))
    fileConn <- file(filename)
  indent <- function(str, indent) {
    return (paste0(paste0(rep(' ', indent), collapse = ''), str))
  }
  
  lines <- c()
  Y_def <- paste0("Y: ", obj$outcome_type)
  S_def <- paste0("S: ", paste0(strtoi(obj$strata, base = 2), collapse = ' '))
  lines <- c(lines, Y_def, S_def)
  lines <- c(lines, "parameter {")
  for (i in 1:length(obj$parameter_list)) {
    if(obj$parameter_list[[i]]$type == "real_vct") {
      type_str <- paste0('vector[', obj$parameter_list[[i]]$dim, ']')
    }
    else
      type_str <- obj$parameter_list[[i]]$type
    name_str <- paste0(
      c(obj$parameter_list[[i]]$class,
      paste0(
        obj$parameter_list[[i]]$group,
        collapse = '_'
      ),
      obj$parameter_list[[i]]$name),
      collapse = '_'
    )
    prior <- obj$priors[[obj$parameter_list[[i]]$prior_type]]
    call <- prior$call
    if (is.null(call))
      tmp_str <- "uniform()"
    else
      tmp_str <- deparse(call, width.cutoff = 500L)
    lines <- c(
      lines, 
      indent(paste0('<', type_str, '> ', name_str, ' ~ ', tmp_str), indent = 4)
    )
  }
  lines <- c(lines, "}")
  
  # Strata models
  lines <- c(lines, "strata {")
  lines <- c(lines, indent(
    paste0(
      strtoi(obj$strata[1], 2), ": 0"
    ), indent = 4
  ))
  for (stratum in strata[-1]){
    prse_fml_S <- parse.formula(obj$S.model)
    name_intrcpt <- paste('S', stratum, "Intrcpt", sep = '_')
    name_coef <- paste('S', stratum, "Coef", sep = '_')
    str_core <- c()
    if (prse_fml_S$has_intercept)
      str_core <- c(str_core, name_intrcpt)
    str_core <- c(str_core, paste0("$ * ", name_coef))
    str_mean <- paste(str_core, collapse = ' + ')
    lines <- c(
      lines, 
      indent(paste0(strtoi(stratum, 2), ": ", str_mean), indent = 4)
    )
  }
  lines <- c(lines, "}")
  # Outcome models
  lines <- c(lines, "outcome {")
  for (stratum in strata){
    for (j in 0:1){
      str_index <- paste0(strtoi(stratum, 2), ", ", j, ": ")
      if (stratum %in% ER)
        group <- c(stratum, 0)
      else 
        group <- c(stratum, j)
      str_model <- get_dist_str(obj$Y.family, obj$Y.model, group)
      lines <- c(
        lines, 
        indent(paste0(str_index, str_model), indent = 4)
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

write.pso(obj, "survival.pso")
