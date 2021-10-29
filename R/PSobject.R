PSObject <- function(
  S.formula, Y.formula, Y.family, data, monotonicity, ER, trunc,
  prior_intercept = prior_uniform(),
  prior_coefficient = prior_normal(),
  prior_sigma = prior_inv_gamma(),
  prior_alpha = prior_inv_gamma(),
  prior_lambda = prior_inv_gamma(),
  prior_theta = prior_normal()
) {
  strata <- get_strata(monotonicity)
  parameter_list <- list()
  
  # models for stratum
  prse_fml_S <- parse.formula(S.formula, data)
  print("haha1")
  
  for (stratum in strata[-1]) {
    if (prse_fml_S$has_intercept) # S has intercept
      parameter_list <- append(parameter_list, list(list(
        type = "real", prior_type = "prior_intercept", name = "Intrcpt",
        group = stratum, class = 'S'
      )))
    parameter_list <- append(parameter_list, list(list(
      type = "real_vct", dim = prse_fml_S$num_of_predictors,
      prior_type = "prior_coefficient", name = "Coef",
      group = stratum, class = 'S'
    )))
  }
  
  # models for outcome
  prse_fml_Y <- parse.formula(Y.formula, data)
  print("haha2")
  Y_params <- get_param_from_model(model(Y.family, Y.formula), data)
  
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
  
  df <- list(
    N = nrow(data), 
    Z = dplyr::pull(data, prse_fml_S$response[1]),
    D = dplyr::pull(data, prse_fml_S$response[2]),
    Y = dplyr::pull(data, prse_fml_Y$response[1]),
    XS = prse_fml_S$model_matrix,
    XY = prse_fml_Y$model_matrix
  )
  if (!is.na(prse_fml_Y$response[2]))
    df$C <- dplyr::pull(data, prse_fml_S$response[2])
  
  return (structure(list(
    PSsettings = list(
      S.formula = S.formula, 
      Y.formula = Y.formula, 
      Y.family = Y.family,
      Y.type = get_outcome_type(Y.family),
      strata = strata,
      ER = ER
    ),
    PSvars = list(
      treatment = prse_fml_S$response[1],
      intervention = prse_fml_S$response[2],
      outcome = prse_fml_Y$response[1],
      censor = prse_fml_Y$response[2]
    ),
    PScovs = list(
      S.intercept = prse_fml_S$has_intercept,
      S.predictors = prse_fml_S$predictors,
      S.dim = prse_fml_S$num_of_predictors,
      Y.intercept = prse_fml_Y$has_intercept,
      Y.predictors = prse_fml_Y$predictors,
      Y.dim = prse_fml_Y$num_of_predictors,
      ncount = nrow(data)
    ),
    data = df,
    parameter_list = parameter_list,
    priors = list(
      prior_intercept = prior_intercept,
      prior_coefficient = prior_coefficient,
      prior_sigma = prior_sigma,
      prior_alpha = prior_alpha,
      prior_lambda = prior_lambda,
      prior_theta = prior_theta
    )
  ), class = "PSObject"))
}

write.pso <- function(obj, filename = NULL){
  if (!is.null(filename))
    fileConn <- file(filename)
  indent <- function(str, indent) {
    return (paste0(paste0(rep(' ', indent), collapse = ''), str))
  }
  
  lines <- c()
  Y_def <- paste0("Y: ", obj$PSsettings$Y.type)
  S_def <- paste0("S: ", paste0(
    strtoi(obj$PSsettings$strata, base = 2), collapse = ' '))
  Dim_def <- paste0("Dim: ", obj$PScovs$S.dim, " ", obj$PScovs$Y.dim)
  lines <- c(lines, Y_def, S_def, Dim_def)
  lines <- c(lines, "parameter {")
  for (i in 1:length(obj$parameter_list)) {
    if(obj$parameter_list[[i]]$type == "real_vct" && obj$parameter_list[[i]]$dim == 0)
      next
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
      strtoi(obj$PSsettings$strata[1], 2), ": 0"
    ), indent = 4
  ))
  for (stratum in obj$PSsettings$strata[-1]){
    name_intrcpt <- paste('S', stratum, "Intrcpt", sep = '_')
    name_coef <- paste('S', stratum, "Coef", sep = '_')
    str_core <- c()
    if (obj$PScovs$S.intercept)
      str_core <- c(str_core, name_intrcpt)
    if (obj$PScovs$S.dim)
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
  for (stratum in obj$PSsettings$strata){
    for (j in 0:1){
      str_index <- paste0(strtoi(stratum, 2), ", ", j, ": ")
      if (stratum %in% obj$PSsettings$ER)
        group <- c(stratum, 0)
      else 
        group <- c(stratum, j)
      str_model <- get_dist_str(obj, group)
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
