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
  else {
    mean_surv_prob <- rstan::extract(PStrataObj$post_samples)$'mean_surv_prob'
    S_count <- PStrataObj$PSobject$strata_info$num_strata
    Z_count <- PStrataObj$PSobject$strata_info$num_treatment
    time_points <- make_standata(PStrataObj$PSobject)$time
    T_count <- length(time_points)
    iter_count <- dim(mean_surv_prob)[1]
    SZDG_table <- PStrataObj$PSobject$SZDG_table
    
    outcome_array <- array(NA, dim = c(S_count, Z_count, T_count, iter_count))
    for (S in 1:S_count) {
      for (Z in 1:Z_count) {
        G <- SZDG_table[SZDG_table[, 'S'] == S - 1 & SZDG_table[, 'Z'] == Z - 1, 'G']
        outcome_array[S, Z, , ] <- t(mean_surv_prob[, G, ])
      }
    }
    
    
    dimnames(outcome_array)[[1]] <- PStrataObj$PSobject$strata_info$strata_names
    dimnames(outcome_array)[[2]] <- PStrataObj$PSobject$Z_names
    dimnames(outcome_array)[[3]] <- 1:T_count
    
    return (structure(
      list(
        outcome_array = outcome_array,
        is.survival = T,
        time_points = time_points
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
  } else {
    cat("Survival PSOutcome Object (", dim(PSoutcome$outcome_array)[1],
        " strata, ", dim(PSoutcome$outcome_array)[2], 
        " treatment arms, ", dim(PSoutcome$outcome_array)[4], " iterations)\n", sep = '')
    cat("Evaluated at ", dim(PSoutcome$outcome_array)[3], " time points.\n", sep = '')
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
    SZ <- purrr::transpose(stringr::str_split(rownames(summary_res_df), ":"))
    summary_res_df <- cbind(S = unlist(SZ[[1]]), Z = unlist(SZ[[2]]), summary_res_df)
    
    if (match.arg(type, c("array", "matrix", "data.frame")) == "array")
      return (summary_res_array)
    else if (match.arg(type, c("array", "matrix", "data.frame")) == "matrix")
      return (summary_res_matrix)
    else if (match.arg(type, c("array", "matrix", "data.frame")) == "data.frame")
      return (summary_res_df)
  }
  else {
    summary_raw <- apply(PSoutcome$outcome_array, c(3,2,1), 
                         function(x) c(mean(x), sd(x), quantile(x, 0.025), quantile(x,0.25),
                                       median(x), quantile(x, 0.75), quantile(x, 0.975)))
    summary_res_array <- aperm(summary_raw, c(4,3,2,1))
    summary_names <- c("mean", "sd", "2.5%", "25%", "median", "75%", "97.5%")
    S_names <- dimnames(summary_res_array)[[1]]
    Z_names <- dimnames(summary_res_array)[[2]]
    T_names <- dimnames(summary_res_array)[[3]]
    dimnames(summary_res_array)[[4]] <- summary_names
    summary_res_matrix <- matrix(nrow = length(S_names) * length(Z_names) * length(T_names), 
                                 ncol = length(summary_names))
    for (i in 1:dim(summary_res_array)[1])
      for (j in 1:dim(summary_res_array)[2]) {
        for (t in 1:dim(summary_res_array)[3]) {
          summary_res_matrix[(i - 1) * dim(summary_res_array)[2] * dim(summary_res_array)[3] + (j - 1) * dim(summary_res_array)[3] + t, ] <- summary_res_array[i, j, t, ]
        }
      }
    tmp_names <- expand.grid(T_names, Z_names, S_names)[, c(3,2,1)]
    colnames(summary_res_matrix) <- summary_names
    rownames(summary_res_matrix) <- apply(tmp_names, 1, paste, collapse = ":")
    
    summary_res_df <- as.data.frame(summary_res_matrix)
    SZT <- purrr::transpose(stringr::str_split(rownames(summary_res_df), ":"))
    if (any(is.na(as.integer(unlist(SZT[[3]]))))){
      summary_res_df <- cbind(S = unlist(SZT[[1]]), Z = unlist(SZT[[2]]), 
                       T = unlist(SZT[[3]]), summary_res_df)
    }
    else{
      summary_res_df <- cbind(S = unlist(SZT[[1]]), Z = unlist(SZT[[2]]), 
                       T = PSoutcome$time_points[as.integer(unlist(SZT[[3]]))], summary_res_df)
    }
    
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
  } else {
    PSoutcome$outcome_array <- PSoutcome$outcome_array[S, Z, T, iter, drop = FALSE]
  }
  
  return (PSoutcome)
}

#' @export
plot.PSOutcome <- function(PSoutcome, se = T) {
  if (!PSoutcome$is.survival) {
    plot_df <- summary(PSoutcome, "data.frame")
    Gplot <- ggplot2::ggplot(plot_df) + ggplot2::geom_point(ggplot2::aes(x = mean, y = "")) + 
      ggplot2::facet_grid(Z~S, scale = "free")
    if (se)
      Gplot <- Gplot + ggplot2::geom_linerange(ggplot2::aes(xmin = `2.5%`, xmax = `97.5%`, y = ""))
    return (Gplot)
  }
  else {
    plot_df <- summary(PSoutcome, "data.frame")
    if (any(is.na(as.integer(plot_df$T)))) {
      Gplot <- ggplot2::ggplot(plot_df) + ggplot2::geom_point(ggplot2::aes(x = T, y = mean)) + 
        ggplot2::facet_grid(Z~S, scale = "free")
      if (se)
        Gplot <- Gplot + ggplot2::geom_linerange(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`, x = T))
      Gplot <- Gplot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
    }
    else{
      Gplot <- ggplot2::ggplot(plot_df) + ggplot2::geom_line(ggplot2::aes(x = T, y = mean)) + 
        ggplot2::facet_grid(Z~S, scale = "free")
      if (se)
        Gplot <- Gplot + ggplot2::geom_ribbon(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`, x = T),
                                              alpha = 0.3)
    }
    return (Gplot)
  }
}
