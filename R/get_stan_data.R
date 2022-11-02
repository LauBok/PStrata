get.stan.data <- function(S.formula, Y.formula, data){
  
  prse_fml_S <- parse.formula(S.formula, data)
  prse_fml_Y <- parse.formula(Y.formula, data)
  length_S <- length(prse_fml_S$response)

  D_bin <- dplyr::pull(data, prse_fml_S$response[2])
  if(length_S > 2){
    for(i in 3:length_S){
      D_bin <- 2 * D_bin + dplyr::pull(data, prse_fml_S$response[i])
    }
  }

  df <- list(
    N = nrow(data), 
    PS = ncol(prse_fml_S$model_matrix),
    PG = ncol(prse_fml_Y$model_matrix),
    Z = dplyr::pull(data, prse_fml_S$response[1]),
    D = D_bin,
    Y = dplyr::pull(data, prse_fml_Y$response[1]),
    XS = prse_fml_S$model_matrix,
    XG = prse_fml_Y$model_matrix
  )
  # censoring
  if (length(prse_fml_Y$response) == 2) {
    df$delta = dplyr::pull(data, prse_fml_Y$response[2])
  }
  
  S_re_list <- prse_fml_S$random_matrix_list
  if (!is.null(S_re_list)) {
    for (i in 1:length(S_re_list)) {
      df[[paste0("PS_RE_", i)]] <- length(S_re_list[[i]]$terms)
      df[[paste0("NS_RE_", i)]] <- length(attr(S_re_list[[i]]$factors, "levels"))
      df[[paste0("XS_RE_", i)]] <- S_re_list[[i]]$matrix
    }
  }
  
  Y_re_list <- prse_fml_Y$random_matrix_list
  if (!is.null(Y_re_list)) {
    for (i in 1:length(Y_re_list)) {
      df[[paste0("PG_RE_", i)]] <- length(Y_re_list[[i]]$terms)
      df[[paste0("NG_RE_", i)]] <- length(attr(Y_re_list[[i]]$factors, "levels"))
      df[[paste0("XG_RE_", i)]] <- Y_re_list[[i]]$matrix
    }
  }
  
  return(df)
}
