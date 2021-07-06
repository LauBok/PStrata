write.pso <- function(obj, filename){
  fileConn <- file(filename)
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
      lines <- c(lines, indent(paste0('$', i, ': ', obj$symbol_list[i]), indent = 4))
    }
  }
  lines <- c(lines, "}")
  # Parameters
  lines <- c(lines, "parameter {")
  for (i in 1:length(obj$param_list)){
    if (obj$param_list[[i]]$model == "normal"){
      type <- "continuous"    
    }
    else if (obj$param_list[[i]]$model == "inv_gamma"){
      type <- "positive"
    }
    lines <- c(lines, indent(
      paste0('@', i, ': <', type, '> ', obj$param_list[[i]]$name),
      indent = 4))
  }
  lines <- c(lines, "}")
  # Prior Distribution
  lines <- c(lines, "prior {")
  for (i in 1:length(obj$param_list)){
    lines <- c(lines, indent(
      paste0('@', i, ' ~ ', obj$param_list[[i]]$model, '(', 
             paste0(obj$param_list[[i]]$args, collapse = ', '), ')'),
      indent = 4
    ))
  }
  lines <- c(lines, "}")
  # Strata models
  lines <- c(lines, "strata {")
  for (i in 1:length(obj$S_list)){
    lines <- c(lines, indent(
      paste0(strtoi(names(obj$S_list)[i], 2), ": ", obj$S_list[[i]]$name),
      indent = 4
    ))
  }
  lines <- c(lines, "}")
  # Outcome models
  lines <- c(lines, "outcome {")
  for (i in 1:length(obj$Y_list)){
    for (j in 0:1){
      lines <- c(lines, indent(
        paste0(strtoi(names(obj$Y_list)[i], 2), ", ", j,  ": ", obj$Y_list[[i]][[j + 1]]$name),
        indent = 4
      ))
    }
  }
  lines <- c(lines, "}")
  writeLines(lines, fileConn)
  close(fileConn)
}
