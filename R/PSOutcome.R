#' Estimated potential outcome for principal stratification analysis
#' 
#' Create an object useful to present the potential outcomes under each treatment arm
#' for each principal stratum. Contrasts between treatment arms or principal strata
#' are easy to obtain from this object.
#' 
#' @param PStrataObj an object of class \code{PStrata} or \code{PStrata_survival}
#' 
#' @export
PSOutcome <- function(PStrataObj) {
  if (class(PStrataObj) == "PStrata"){
    mean_effect <- rstan::extract(PStrataObj$post_samples)$'mean_effect'
    S_count <- PStrataObj$PSobject$strata_info$num_strata
    Z_count <- PStrataObj$PSobject$strata_info$num_treatment
    iter_count <- nrow(mean_effect)
    SZDG_table <- PStrataObj$PSobject$SZDG_table
    
    outcome_array <- array(NA, dim = c(S_count, Z_count, iter_count))
    for (S in 1:S_count) {
      for (Z in 1:Z_count) {
        G <- SZDG_table[SZDG_table[, 'S'] == S - 1 & SZDG_table[, 'Z'] == Z - 1, 'G']
        outcome_array[S, Z, ] <- mean_effect[, G]
      }
    }
    
    dimnames(outcome_array)[[1]] <- PStrataObj$PSobject$strata_info$strata_names
    dimnames(outcome_array)[[2]] <- PStrataObj$PSobject$Z_names
    
    return (structure(
      list(
        outcome_array = outcome_array,
        is.survival = F
      ),
      class = "PSOutcome"
    ))
  }
    
}

#' @exportS3Method 
print.PSOutcome <- function(PSoutcome) {
  if (!PSoutcome$is.survival) {
    cat("Non-survival PSOutcome Object (", dim(PSoutcome$outcome_array)[1],
        " strata, ", dim(PSoutcome$outcome_array)[2], 
        " treatment arms, ", dim(PSoutcome$outcome_array)[3], " iterations)\n", sep = '')
  }
}

#' @export
summary.PSOutcome <- function(PSoutcome, type = c("array", "matrix", "data.frame")){
  if (!PSoutcome$is.survival) {
    summary_raw <- apply(PSoutcome$outcome_array, c(2,1), 
                         function(x) c(mean(x), sd(x), quantile(x, 0.025), quantile(x,0.25),
                                       median(x), quantile(x, 0.75), quantile(x, 0.975)))
    summary_res_array <- aperm(summary_raw, c(3,2,1))
    summary_names <- c("mean", "sd", "2.5%", "25%", "median", "75%", "97.5%")
    S_names <- dimnames(summary_res_array)[[1]]
    Z_names <- dimnames(summary_res_array)[[2]]
    dimnames(summary_res_array)[[3]] <- summary_names
    summary_res_matrix <- matrix(nrow = length(S_names) * length(Z_names), 
                                 ncol = length(summary_names))
    for (i in 1:dim(summary_res_array)[1])
      for (j in 1:dim(summary_res_array)[2]) {
        summary_res_matrix[(i - 1) * dim(summary_res_array)[2] + j, ] <- summary_res_array[i, j, ]
      }
    tmp_names <- expand.grid(Z_names, S_names)[, c(2,1)]
    colnames(summary_res_matrix) <- summary_names
    rownames(summary_res_matrix) <- apply(tmp_names, 1, paste, collapse = ":")
    
    summary_res_df <- as.data.frame(summary_res_matrix)
    
    if (match.arg(type, c("array", "matrix", "data.frame")) == "array")
      return (summary_res_array)
    else if (match.arg(type, c("array", "matrix", "data.frame")) == "matrix")
      return (summary_res_matrix)
    else if (match.arg(type, c("array", "matrix", "data.frame")) == "data.frame")
      return (summary_res_df)
  }
  
}

#' @export
`[.PSOutcome` <- function(PSoutcome, S, Z, T, iter) {
  if (missing(S)) S <- TRUE
  if (missing(Z)) Z <- TRUE
  if (missing(T)) T <- TRUE
  if (missing(iter)) iter <- TRUE
  if (!PSoutcome$is.survival){
    if (!all(iter))
      PSoutcome$outcome_array <- PSoutcome$outcome_array[S, Z, iter, drop = FALSE]
    else
      PSoutcome$outcome_array <- PSoutcome$outcome_array[S, Z, T, drop = FALSE]
  }
  
  return (PSoutcome)
}

#' @export
plot.PSOutcome <- function(PSoutcome, se = T) {
  if (!PSoutcome$is.survival) {
    plot_df <- summary(PSoutcome, "data.frame")
    SZ <- purrr::transpose(stringr::str_split(rownames(plot_df), ":"))
    plot_df <- cbind(S = unlist(SZ[[1]]), Z = unlist(SZ[[2]]), plot_df)
    Gplot <- ggplot2::ggplot(plot_df) + ggplot2::geom_point(ggplot2::aes(x = mean, y = "")) + 
      ggplot2::facet_grid(Z~S, scale = "free")
    if (se)
      Gplot <- Gplot + ggplot2::geom_linerange(ggplot2::aes(xmin = `2.5%`, xmax = `97.5%`, y = ""))
    return (Gplot)
  }
}
