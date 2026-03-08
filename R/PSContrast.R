#' Contrast of potential outcome for principal stratification analysis
#'
#' Create an object that represents contrast of potential outcomes by treatment arms,
#' strata or time points.
#'
#' @param outcome an object of class \code{PSoutcome} or \code{PSContrast}
#' @param S a vector denoting which strata to take contrasts. Default is \code{NULL} indicating
#' no contrasts are taken. Set to `TRUE` to take contrasts between all strata.
#' @param Z a vector denoting which treatment arms to take contrasts. Default is \code{NULL} indicating
#' no contrasts are taken. Set to `TRUE` to take contrasts between all treatment arms.
#' @param T a vector denoting which time points to take contrasts. Default is \code{NULL} indicating
#' no contrasts are taken. Set to `TRUE` to take contrasts between all time points. This is used
#' only when `object` is obtained under survival outcome.
#' @param type Either \code{"all"} (default), \code{"sequential"} or \code{"cycle"}.
#' If \code{"all"}, every pairwise contrasts are taken.
#' If \code{"sequential"}, contrasts are taken over every consecutive pairs.
#' If \code{"cycle"}, contrasts are taken over every consecutive pairs and also between the
#' first and the last levels.
#'
#' @return An S3 object of class \code{PSContrast} and \code{PSOutcome}, containing
#' \item{outcome_array}{A num_strata * num_treatment * num_iter array of contrast if the outcome type is non-survival
#' or a num_strata * num_treatment * num_time_points * num_iter array of contrast if the outcome type is survival.}
#' \item{is.survival}{A boolean value, whether the outcome type is survival.}
#' \item{time_points}{The time points at which the outcome is evaluated, if the outcome type is survival.}
#' The S3 method \code{summary} and \code{plot} can be applied to the returned object.
#'
#' @export
PSContrast <- function(outcome, S = NULL, Z = NULL, T = NULL,
                       type = c("all", "sequential", "cycle")) {
  type <- match.arg(type)
  arr <- outcome$outcome_array

  # Apply contrasts along each requested dimension
  if (!is.null(S)) arr <- contrast_along_dim(arr, 1L, S, type)
  if (!is.null(Z)) arr <- contrast_along_dim(arr, 2L, Z, type)
  if (!is.null(T) && outcome$is.survival) arr <- contrast_along_dim(arr, 3L, T, type)

  structure(
    list(
      outcome_array = arr,
      is.survival   = outcome$is.survival,
      time_points   = outcome$time_points
    ),
    class = c("PSContrast", "PSOutcome")
  )
}

#' @export
print.PSContrast <- function(x, ...) {
  dims <- dim(x$outcome_array)
  if (!x$is.survival) {
    cat("Non-survival PSContrast Object (", dims[1],
        " strata, ", dims[2], " treatment arms, ", dims[3], " iterations)\n", sep = "")
  } else {
    cat("Survival PSContrast Object (", dims[1],
        " strata, ", dims[2], " treatment arms, ", dims[4], " iterations)\n", sep = "")
    cat("Evaluated at ", dims[3], " time points.\n", sep = "")
  }
}
