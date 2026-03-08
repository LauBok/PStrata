#' Estimated potential outcome for principal stratification analysis
#'
#' Create an object useful to present the potential outcomes under each treatment arm
#' for each principal stratum. Contrasts between treatment arms or principal strata
#' are easy to obtain from this object.
#'
#' @param PStrataObj an object of class \code{PStrata} or \code{PStrata_survival}
#' @param type whether the causal estimand is survival probability or RACE, ignored for non-survival outcomes.
#' @return An S3 object of type \code{PSOutcome}, containing
#' \item{outcome_array}{A num_strata * num_treatment * num_iter array of mean outcome if the outcome type is non-survival
#' or a num_strata * num_treatment * num_time_points * num_iter array of mean outcome if the outcome type is survival.}
#' \item{is.survival}{A boolean value, whether the outcome type is survival.}
#' \item{time_points}{The time points at which the outcome is evaluated, if the outcome type is survival.}
#' The S3 method \code{summary} and \code{plot} can be applied to the returned object.
#' @export
PSOutcome <- function(PStrataObj, type = c("probability", "RACE")) {
  info <- PStrataObj$PSobject
  SZDG_table <- info$SZDG_table
  S_count <- info$strata_info$num_strata
  Z_count <- info$strata_info$num_treatment

  if (inherits(PStrataObj, "PStrata") && !inherits(PStrataObj, "PStrata_survival")) {
    outcome_array <- build_outcome_array_ordinary(PStrataObj, S_count, Z_count, SZDG_table)
    dimnames(outcome_array)[[1]] <- info$strata_info$strata_names
    dimnames(outcome_array)[[2]] <- info$Z_names

    return(new_PSOutcome(outcome_array, is.survival = FALSE))
  }

  # Survival case
  outcome_type <- match.arg(type, c("probability", "RACE"))
  outcome_array <- build_outcome_array_survival(
    PStrataObj, outcome_type, S_count, Z_count, SZDG_table
  )
  time_points <- make_standata(info)$time

  dimnames(outcome_array)[[1]] <- info$strata_info$strata_names
  dimnames(outcome_array)[[2]] <- info$Z_names
  dimnames(outcome_array)[[3]] <- seq_len(dim(outcome_array)[3])

  new_PSOutcome(outcome_array, is.survival = TRUE, time_points = time_points)
}

# --- Internal constructors ---

#' @noRd
new_PSOutcome <- function(outcome_array, is.survival, time_points = NULL) {
  structure(
    list(
      outcome_array = outcome_array,
      is.survival   = is.survival,
      time_points   = time_points
    ),
    class = "PSOutcome"
  )
}

#' Build the [S, Z, iter] outcome array for non-survival
#' @noRd
build_outcome_array_ordinary <- function(PStrataObj, S_count, Z_count, SZDG_table) {
  mean_effect <- rstan::extract(PStrataObj$post_samples)$'mean_effect'
  iter_count <- nrow(mean_effect)
  arr <- array(NA, dim = c(S_count, Z_count, iter_count))
  for (s in seq_len(S_count)) {
    for (z in seq_len(Z_count)) {
      G <- SZDG_table[SZDG_table[, "S"] == s - 1 & SZDG_table[, "Z"] == z - 1, "G"]
      arr[s, z, ] <- mean_effect[, G]
    }
  }
  arr
}

#' Build the [S, Z, T, iter] outcome array for survival
#' @noRd
build_outcome_array_survival <- function(PStrataObj, outcome_type, S_count, Z_count, SZDG_table) {
  extract_name <- if (outcome_type == "probability") "mean_surv_prob" else "mean_RACE"
  mean_outcome <- rstan::extract(PStrataObj$post_samples)[[extract_name]]
  T_count <- dim(mean_outcome)[3]
  iter_count <- dim(mean_outcome)[1]

  arr <- array(NA, dim = c(S_count, Z_count, T_count, iter_count))
  for (s in seq_len(S_count)) {
    for (z in seq_len(Z_count)) {
      G <- SZDG_table[SZDG_table[, "S"] == s - 1 & SZDG_table[, "Z"] == z - 1, "G"]
      arr[s, z, , ] <- t(mean_outcome[, G, ])
    }
  }
  arr
}

# --- S3 methods ---

#' @exportS3Method
print.PSOutcome <- function(x, ...) {
  dims <- dim(x$outcome_array)
  if (!x$is.survival) {
    cat("Non-survival PSOutcome Object (", dims[1],
        " strata, ", dims[2], " treatment arms, ", dims[3], " iterations)\n", sep = "")
  } else {
    cat("Survival PSOutcome Object (", dims[1],
        " strata, ", dims[2], " treatment arms, ", dims[4], " iterations)\n", sep = "")
    cat("Evaluated at ", dims[3], " time points.\n", sep = "")
  }
}

#' @export
summary.PSOutcome <- function(object, type = c("array", "matrix", "data.frame"), ...) {
  type <- match.arg(type)

  summary_array <- summarize_last_dim(object$outcome_array)
  ndim <- length(dim(summary_array))
  dimnames(summary_array)[[ndim]] <- summary_stat_names()

  if (type == "array") return(summary_array)

  mat <- summary_array_to_matrix(summary_array)
  if (type == "matrix") return(mat)

  col_names <- if (object$is.survival) c("S", "Z", "T") else c("S", "Z")
  summary_matrix_to_df(mat, col_names, object$time_points)
}

#' @export
`[.PSOutcome` <- function(PSoutcome, S, Z, T, iter) {
  if (missing(S)) S <- TRUE
  if (missing(Z)) Z <- TRUE
  if (missing(T)) T <- TRUE
  if (missing(iter)) iter <- TRUE

  if (!PSoutcome$is.survival) {
    if (!all(iter))
      PSoutcome$outcome_array <- PSoutcome$outcome_array[S, Z, iter, drop = FALSE]
    else
      PSoutcome$outcome_array <- PSoutcome$outcome_array[S, Z, T, drop = FALSE]
  } else {
    PSoutcome$outcome_array <- PSoutcome$outcome_array[S, Z, T, iter, drop = FALSE]
  }

  PSoutcome
}

#' @export
plot.PSOutcome <- function(x, se = TRUE, ...) {
  lwr <- upr <- NULL
  plot_df <- summary(x, "data.frame")
  plot_df$lwr <- plot_df$`2.5%`
  plot_df$upr <- plot_df$`97.5%`

  if (!x$is.survival) {
    Gplot <- ggplot2::ggplot(plot_df) +
      ggplot2::geom_point(ggplot2::aes(x = mean, y = "")) +
      ggplot2::facet_grid(Z ~ S, scale = "free")
    if (se)
      Gplot <- Gplot +
        ggplot2::geom_linerange(ggplot2::aes(xmin = lwr, xmax = upr, y = ""))
  } else {
    if (any(is.na(suppressWarnings(as.numeric(plot_df$T))))) {
      Gplot <- ggplot2::ggplot(plot_df) +
        ggplot2::geom_point(ggplot2::aes(x = T, y = mean)) +
        ggplot2::facet_grid(Z ~ S, scale = "free")
      if (se)
        Gplot <- Gplot +
          ggplot2::geom_linerange(ggplot2::aes(ymin = lwr, ymax = upr, x = T))
      Gplot <- Gplot +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
    } else {
      Gplot <- ggplot2::ggplot(plot_df) +
        ggplot2::geom_line(ggplot2::aes(x = T, y = mean)) +
        ggplot2::facet_grid(Z ~ S, scale = "free")
      if (se)
        Gplot <- Gplot +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr, x = T), alpha = 0.3)
    }
  }

  Gplot
}
