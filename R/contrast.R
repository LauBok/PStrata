# ===========================================================================
# contrast.R — Compute causal contrasts from PSEstimate
# ===========================================================================

#' Contrast Potential Outcomes
#'
#' Compute pairwise differences of posterior potential outcomes
#' along strata, treatment arms, or time points.
#'
#' @param object A \code{PSEstimate}, \code{PSContrast}, or
#'   \code{PStrataFit} object. If \code{PStrataFit}, \code{estimate()}
#'   is called automatically.
#' @param S Strata to contrast. \code{TRUE} for all pairwise, or a
#'   numeric/logical index vector. \code{NULL} (default) = no contrast.
#' @param Z Treatment arms to contrast. \code{TRUE} for all pairwise.
#' @param T Time points to contrast (survival only).
#' @param type \code{"all"} (all pairwise), \code{"sequential"}
#'   (consecutive pairs), or \code{"cycle"} (consecutive + wrap-around).
#' @param ... Additional arguments passed to \code{estimate()} when
#'   \code{object} is a \code{PStrataFit}.
#' @return A \code{PSContrast} object (inherits from \code{PSEstimate}).
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
#'
#' # Treatment effect (Z contrast) for each stratum
#' ctr_z <- contrast(ps_fit, Z = TRUE)
#' summary(ctr_z)
#' plot(ctr_z)
#'
#' # Stratum contrasts for each treatment arm
#' ctr_s <- contrast(ps_fit, S = TRUE)
#' summary(ctr_s)
#' }
#' @export
contrast <- function(object, ...) UseMethod("contrast")

#' @rdname contrast
#' @export
contrast.PSEstimate <- function(object, S = NULL, Z = NULL, T = NULL,
                                type = c("all", "sequential", "cycle"),
                                ...) {
  type <- match.arg(type)
  arr <- object$outcome_array

  if (!is.null(S)) arr <- contrast_along_dim(arr, 1L, S, type)
  if (!is.null(Z)) arr <- contrast_along_dim(arr, 2L, Z, type)
  if (!is.null(T) && object$is_survival)
    arr <- contrast_along_dim(arr, 3L, T, type)

  structure(list(
    outcome_array = arr,
    is_survival   = object$is_survival,
    time_points   = object$time_points,
    model_info    = object$model_info
  ), class = c("PSContrast", "PSEstimate"))
}

#' @rdname contrast
#' @export
contrast.PStrataFit <- function(object, S = NULL, Z = NULL, T = NULL,
                                type = c("all", "sequential", "cycle"),
                                ...) {
  est <- estimate(object, ...)
  contrast(est, S = S, Z = Z, T = T, type = type)
}

# ===========================================================================
# PSContrast: print
# ===========================================================================

#' @export
print.PSContrast <- function(x, ...) {
  dims <- dim(x$outcome_array)

  if (!x$is_survival) {
    cat("PSContrast (", dims[1], " x ", dims[2], " x ",
        dims[3], " draws)\n", sep = "")
  } else {
    cat("PSContrast [survival] (", dims[1], " x ", dims[2],
        " x ", dims[3], " time points x ", dims[4],
        " draws)\n", sep = "")
  }

  # Show compact summary: mean and 95% CI
  summary_arr <- summarize_last_dim(x$outcome_array)
  ndim <- length(dim(summary_arr))
  dimnames(summary_arr)[[ndim]] <- summary_stat_names()

  mat <- summary_array_to_matrix(summary_arr)
  show_cols <- c("mean", "2.5%", "97.5%")
  cat("\nPosterior summary:\n")
  print(round(mat[, show_cols, drop = FALSE], 4))

  cat("\nUse summary() for full posterior intervals.\n")
  invisible(x)
}

# ===========================================================================
# PSContrast: plot
# ===========================================================================

#' @export
plot.PSContrast <- function(x, ...) {
  arr <- x$outcome_array
  dims <- dim(arr)
  dn <- dimnames(arr)

  if (!x$is_survival) {
    df <- expand.grid(
      iter = seq_len(dims[3]),
      Z    = dn[[2]],
      S    = dn[[1]],
      stringsAsFactors = FALSE
    )
    df$value <- as.vector(aperm(arr, c(3, 2, 1)))
    df$S <- factor(df$S, levels = dn[[1]])

    if (dims[2] == 1) {
      df$label <- paste0(df$S, ": ", dn[[2]])
    } else {
      df$label <- paste0(df$S, ": ", df$Z)
    }
    df$label <- factor(df$label, levels = unique(df$label))

    ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = "steelblue", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          color = "red") +
      ggplot2::facet_wrap(~ label, scales = "free") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y  = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  } else {
    # Survival contrast: ribbon plot with reference line at 0
    plot_df <- summary(x, "data.frame")
    plot_df$T <- as.numeric(plot_df$T)

    ggplot2::ggplot(plot_df,
        ggplot2::aes(x = T, y = mean)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`),
        alpha = 0.2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "red") +
      ggplot2::facet_wrap(~ S) +
      ggplot2::labs(x = "Time", y = "Contrast") +
      ggplot2::theme_minimal()
  }
}
