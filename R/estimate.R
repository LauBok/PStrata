# ===========================================================================
# estimate.R — Extract posterior potential outcomes from PStrataFit
# ===========================================================================

#' Extract Posterior Potential Outcomes
#'
#' Maps posterior samples of outcome-group means back to
#' the (stratum x treatment) grid via the SZDG table.
#'
#' @param object A \code{PStrataFit} object.
#' @param type For survival models: \code{"probability"} (survival probability)
#'   or \code{"RACE"} (restricted average causal effect).
#' @param ... Additional arguments (unused).
#' @return A \code{PSEstimate} object containing
#' \item{outcome_array}{[S, Z, iter] array for non-survival,
#'   or [S, Z, T, iter] for survival.}
#' \item{is_survival}{Logical.}
#' \item{time_points}{Numeric vector (survival only).}
#' \item{model_info}{List of strata names, treatment info,
#'   family, ER flags, outcome groups.}
#' @examples
#' \donttest{
#' data(sim_data_normal)
#' model <- PStrataModel(
#'   S.formula = Z + D ~ 1,
#'   Y.formula = Y ~ 1,
#'   Y.family  = gaussian(),
#'   strata    = c(n = "00", c = "01", a = "11"),
#'   ER        = c("n", "a")
#' )
#' ps_fit <- fit(model, data = sim_data_normal, chains = 2, iter = 500)
#' est <- estimate(ps_fit)
#' summary(est)
#' plot(est)
#' }
#' @export
estimate <- function(object, ...) UseMethod("estimate")

#' @rdname estimate
#' @export
estimate.PStrataFit <- function(object, type = c("probability", "RACE"), ...) {
  m <- object$model
  SZDG <- object$SZDG_table
  samples <- rstan::extract(object$stanfit)

  S_count <- m$num_strata
  Z_count <- m$num_treatment
  Z_names <- paste0(m$treatment_name, "=", 0:(Z_count - 1))

  mi <- list(
    strata_names          = m$strata_names,
    treatment_name        = m$treatment_name,
    Z_names               = Z_names,
    family                = m$outcome_family,
    exclusion_restriction = m$exclusion_restriction,
    outcome_groups        = m$outcome_groups
  )

  if (!m$is_survival) {
    mean_effect <- samples$mean_effect
    arr <- .map_SZDG_to_array(mean_effect, SZDG, S_count, Z_count)
    dimnames(arr) <- list(S = m$strata_names, Z = Z_names, iter = NULL)

    return(structure(list(
      outcome_array = arr,
      is_survival   = FALSE,
      time_points   = NULL,
      model_info    = mi
    ), class = "PSEstimate"))
  }

  # --- survival ---
  type <- match.arg(type)
  extract_name <- if (type == "probability") "mean_surv_prob" else "mean_RACE"
  mean_outcome <- samples[[extract_name]]
  T_count    <- dim(mean_outcome)[3]
  iter_count <- dim(mean_outcome)[1]

  arr <- array(NA, dim = c(S_count, Z_count, T_count, iter_count))
  for (s in seq_len(S_count)) {
    for (z in seq_len(Z_count)) {
      G <- SZDG[SZDG[, "S"] == s - 1 & SZDG[, "Z"] == z - 1, "G"]
      arr[s, z, , ] <- t(mean_outcome[, G, ])
    }
  }
  dimnames(arr) <- list(
    S = m$strata_names, Z = Z_names,
    T = seq_len(T_count), iter = NULL
  )

  structure(list(
    outcome_array = arr,
    is_survival   = TRUE,
    time_points   = object$time_points,
    type          = type,
    model_info    = mi
  ), class = c("PSEstimateSurvival", "PSEstimate"))
}

# --- internal helper ---

#' Map mean_effect[iter, G] to [S, Z, iter] via SZDG
#' @noRd
.map_SZDG_to_array <- function(mean_effect, SZDG, S_count, Z_count) {
  iter_count <- nrow(mean_effect)
  arr <- array(NA, dim = c(S_count, Z_count, iter_count))
  for (s in seq_len(S_count)) {
    for (z in seq_len(Z_count)) {
      G <- SZDG[SZDG[, "S"] == s - 1 & SZDG[, "Z"] == z - 1, "G"]
      arr[s, z, ] <- mean_effect[, G]
    }
  }
  arr
}

# ===========================================================================
# PSEstimate: print
# ===========================================================================

#' @export
print.PSEstimate <- function(x, ...) {
  dims <- dim(x$outcome_array)
  mi <- x$model_info

  if (!x$is_survival) {
    cat("PSEstimate (", dims[1], " strata x ", dims[2],
        " treatments x ", dims[3], " draws)\n", sep = "")
  } else {
    cat("PSEstimate [survival] (", dims[1], " strata x ", dims[2],
        " treatments x ", dims[3], " time points x ", dims[4],
        " draws)\n", sep = "")
  }

  # Posterior mean table (non-survival only)
  if (!x$is_survival) {
    summary_arr <- summarize_last_dim(x$outcome_array)
    mean_mat <- summary_arr[, , 1, drop = TRUE]
    if (is.null(dim(mean_mat)))
      mean_mat <- matrix(mean_mat, nrow = dims[1], ncol = dims[2])
    rownames(mean_mat) <- mi$strata_names
    colnames(mean_mat) <- mi$Z_names

    cat("\nPosterior mean E[Y(z) | S=s]:\n")
    print(round(mean_mat, 4))
  }

  cat("\nUse summary() for posterior intervals.\n")
  cat("Use contrast() to compute causal effects.\n")
  invisible(x)
}

# ===========================================================================
# PSEstimate: summary
# ===========================================================================

#' @export
summary.PSEstimate <- function(object,
                               type = c("array", "matrix", "data.frame"),
                               ...) {
  type <- match.arg(type)

  summary_array <- summarize_last_dim(object$outcome_array)
  ndim <- length(dim(summary_array))
  dimnames(summary_array)[[ndim]] <- summary_stat_names()

  if (type == "array") return(summary_array)

  mat <- summary_array_to_matrix(summary_array)
  if (type == "matrix") return(mat)

  col_names <- if (object$is_survival) c("S", "Z", "T") else c("S", "Z")
  summary_matrix_to_df(mat, col_names, object$time_points)
}

# ===========================================================================
# PSEstimate: plot
# ===========================================================================

#' @export
plot.PSEstimate <- function(x, ...) {
  mi <- x$model_info

  if (!x$is_survival) {
    arr <- x$outcome_array  # [S, Z, iter]
    dims <- dim(arr)
    dn <- dimnames(arr)

    df <- expand.grid(
      iter = seq_len(dims[3]),
      Z    = dn[[2]],
      S    = dn[[1]],
      stringsAsFactors = FALSE
    )
    df$value <- as.vector(aperm(arr, c(3, 2, 1)))
    df$S <- factor(df$S, levels = dn[[1]])

    df$label <- paste0("E[Y|", df$S, ",", df$Z, "]")
    df$label <- factor(df$label, levels = unique(df$label))

    ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = "steelblue", alpha = 0.5) +
      ggplot2::facet_wrap(~ label, scales = "free") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y  = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())

  } else {
    # Survival: ribbon plot over time, faceted by stratum
    plot_df <- summary(x, "data.frame")
    plot_df$T <- as.numeric(plot_df$T)

    ggplot2::ggplot(plot_df,
        ggplot2::aes(x = T, y = mean, color = Z, fill = Z)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`),
        alpha = 0.2, color = NA) +
      ggplot2::facet_wrap(~ S) +
      ggplot2::labs(x = "Time", y = "Survival probability",
                    color = "Treatment", fill = "Treatment") +
      ggplot2::theme_minimal()
  }
}

# ===========================================================================
# PSEstimate: subsetting
# ===========================================================================

#' @export
`[.PSEstimate` <- function(x, S, Z, T, iter) {
  if (missing(S)) S <- TRUE
  if (missing(Z)) Z <- TRUE
  if (missing(T)) T <- TRUE
  if (missing(iter)) iter <- TRUE

  if (!x$is_survival) {
    if (!isTRUE(iter))
      x$outcome_array <- x$outcome_array[S, Z, iter, drop = FALSE]
    else
      x$outcome_array <- x$outcome_array[S, Z, T, drop = FALSE]
  } else {
    x$outcome_array <- x$outcome_array[S, Z, T, iter, drop = FALSE]
  }

  x
}
