#' Data for \pkg{PStrata} Models
#'
#' Generate data for \pkg{PStrata} models to be passed to \bold{Stan}
#' 
#' @param PSobject an object of class \code{PSObject}
#' 
#' @returns a named list of objects containing the required data to fit a
#' \pkg{PStrata} model with \bold{Stan}.
#' 
#' 
#' @export
make_standata <- function(PSobject){
  stopifnot("Rows of the data frames provided for S.formula and Y.formula must macth" =
              nrow(PSobject$S.formula$data) == nrow(PSobject$Y.formula$data))
  
  D_names <- PSobject$S.formula$response_names[-1]
  D_data <- dplyr::select(PSobject$S.formula$data, dplyr::all_of(c("D")))
  D_int <- apply(D_data, 1, 
        function(x) 
          t(c(x, 0)) %*% c(1, cumprod(PSobject$strata_info$max_postrand_level + 1)))
  
  Z <- dplyr::pull(PSobject$S.formula$data, PSobject$S.formula$response_names[1])
  if (is.factor(Z)) {
    Z_int <- as.integer(Z) - 1
  } else if (all(sapply(Z, function(x) x %in% 0:(PSobject$strata_info$num_treatment - 1)))) {
    Z_int <- Z
  } else {
    Z_int <- as.integer(as.factor(Z)) - 1
  }
  
  df <- list(
    N = nrow(PSobject$S.formula$data), 
    PS = ncol(PSobject$S.formula$fixed_eff_matrix),
    PG = ncol(PSobject$Y.formula$fixed_eff_matrix),
    Z = Z_int,
    D = D_int,
    Y = dplyr::pull(PSobject$Y.formula$data, PSobject$Y.formula$response_names[1]),
    XS = PSobject$S.formula$fixed_eff_matrix,
    XG = PSobject$Y.formula$fixed_eff_matrix
  )
  # censoring
  if (PSobject$is.survival) {
    df$delta = dplyr::pull(PSobject$Y.formula$data, PSobject$Y.formula$response_names[2])
    if (length(PSobject$survival.time.points) == 1) {
      df$time <- seq(0, max(df$Y) * 0.9, length.out = PSobject$survival.time.points)
    }
    else{
      df$time <- PSobject$survival.time.points
    }
    df$T <- length(df$time)
  }
  
  S_re_list <- PSobject$S.formula$random_eff_list
  if (!is.null(S_re_list)) {
    for (i in 1:length(S_re_list)) {
      df[[paste0("PS_RE_", i)]] <- length(S_re_list[[i]]$terms)
      df[[paste0("NS_RE_", i)]] <- length(attr(S_re_list[[i]]$factors, "levels"))
      df[[paste0("XS_RE_", i)]] <- S_re_list[[i]]$matrix
    }
  }
  
  Y_re_list <- PSobject$Y.formula$random_eff_list
  if (!is.null(Y_re_list)) {
    for (i in 1:length(Y_re_list)) {
      df[[paste0("PG_RE_", i)]] <- length(Y_re_list[[i]]$terms)
      df[[paste0("NG_RE_", i)]] <- length(attr(Y_re_list[[i]]$factors, "levels"))
      df[[paste0("XG_RE_", i)]] <- Y_re_list[[i]]$matrix
    }
  }
  
  return(df)
}
