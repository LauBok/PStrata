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
  old_outcome_array <- outcome$outcome_array
  type <- match.arg(type, c("all", "sequential", "cycle"))
  if (outcome$is.survival == F) {
    if (!is.null(S)) {
      dim_outcome_array <- dim(old_outcome_array)
      name_outcome_array <- dimnames(old_outcome_array)
      levels <- (1:dim_outcome_array[1])[S]
      if (type == "all")
        contrast_levels <- (function(x) x[x[,1] > x[,2], ])(expand.grid(levels, levels))
      if (type == "sequential" || (type == "cycle" && length(levels) == 2))
        contrast_levels <- data.frame(Var1 = levels[-1], Var2 = levels[-length(levels)])
      if (type == "cycle" && length(levels) > 2) {
        contrast_levels <- data.frame(Var1 = c(levels[-1], levels[1]), Var2 = levels)
      }
      dim_outcome_array[1] <- nrow(contrast_levels)
      name_outcome_array[[1]] <- apply(contrast_levels, 1, function(x) {
        paste0("{", name_outcome_array[[1]][x[1]], "}-{", name_outcome_array[[1]][x[2]], "}")
      })
      new_outcome_array <- array(NA, dim = dim_outcome_array, dimnames = name_outcome_array)
      for (i in 1:nrow(contrast_levels)) {
        new_outcome_array[i, , ] <- old_outcome_array[contrast_levels[i, 1], ,] - 
          old_outcome_array[contrast_levels[i, 2], ,]
      }
      old_outcome_array <- new_outcome_array
    }
    if (!is.null(Z)) {
      dim_outcome_array <- dim(old_outcome_array)
      name_outcome_array <- dimnames(old_outcome_array)
      levels <- (1:dim_outcome_array[2])[Z]
      if (type == "all")
        contrast_levels <- (function(x) x[x[,1] > x[,2], ])(expand.grid(levels, levels))
      if (type == "sequential" || (type == "cycle" && length(levels) == 2))
        contrast_levels <- data.frame(Var1 = levels[-1], Var2 = levels[-length(levels)])
      if (type == "cycle" && length(levels) > 2) {
        contrast_levels <- data.frame(Var1 = c(levels[-1], levels[1]), Var2 = levels)
      }
      dim_outcome_array[2] <- nrow(contrast_levels)
      name_outcome_array[[2]] <- apply(contrast_levels, 1, function(x) {
        paste0("{", name_outcome_array[[2]][x[1]], "}-{", name_outcome_array[[2]][x[2]], "}")
      })
      new_outcome_array <- array(NA, dim = dim_outcome_array, dimnames = name_outcome_array)
      for (i in 1:nrow(contrast_levels)) {
        new_outcome_array[, i, ] <- old_outcome_array[, contrast_levels[i, 1] ,] - 
          old_outcome_array[, contrast_levels[i, 2], ]
      }
      old_outcome_array <- new_outcome_array
    }
  } else {
    if (!is.null(S)) {
      dim_outcome_array <- dim(old_outcome_array)
      name_outcome_array <- dimnames(old_outcome_array)
      levels <- (1:dim_outcome_array[1])[S]
      if (type == "all")
        contrast_levels <- (function(x) x[x[,1] > x[,2], ])(expand.grid(levels, levels))
      if (type == "sequential" || (type == "cycle" && length(levels) == 2))
        contrast_levels <- data.frame(Var1 = levels[-1], Var2 = levels[-length(levels)])
      if (type == "cycle" && length(levels) > 2) {
        contrast_levels <- data.frame(Var1 = c(levels[-1], levels[1]), Var2 = levels)
      }
      dim_outcome_array[1] <- nrow(contrast_levels)
      name_outcome_array[[1]] <- apply(contrast_levels, 1, function(x) {
        paste0("{", name_outcome_array[[1]][x[1]], "}-{", name_outcome_array[[1]][x[2]], "}")
      })
      new_outcome_array <- array(NA, dim = dim_outcome_array, dimnames = name_outcome_array)
      for (i in 1:nrow(contrast_levels)) {
        new_outcome_array[i, , ,] <- old_outcome_array[contrast_levels[i, 1], , ,] - 
          old_outcome_array[contrast_levels[i, 2], , ,]
      }
      old_outcome_array <- new_outcome_array
    }
    if (!is.null(Z)) {
      dim_outcome_array <- dim(old_outcome_array)
      name_outcome_array <- dimnames(old_outcome_array)
      levels <- (1:dim_outcome_array[2])[Z]
      if (type == "all")
        contrast_levels <- (function(x) x[x[,1] > x[,2], ])(expand.grid(levels, levels))
      if (type == "sequential" || (type == "cycle" && length(levels) == 2))
        contrast_levels <- data.frame(Var1 = levels[-1], Var2 = levels[-length(levels)])
      if (type == "cycle" && length(levels) > 2) {
        contrast_levels <- data.frame(Var1 = c(levels[-1], levels[1]), Var2 = levels)
      }
      dim_outcome_array[2] <- nrow(contrast_levels)
      name_outcome_array[[2]] <- apply(contrast_levels, 1, function(x) {
        paste0("{", name_outcome_array[[2]][x[1]], "}-{", name_outcome_array[[2]][x[2]], "}")
      })
      new_outcome_array <- array(NA, dim = dim_outcome_array, dimnames = name_outcome_array)
      for (i in 1:nrow(contrast_levels)) {
        new_outcome_array[, i, , ] <- old_outcome_array[, contrast_levels[i, 1] , ,] - 
          old_outcome_array[, contrast_levels[i, 2], , ]
      }
      old_outcome_array <- new_outcome_array
    }
    if (!is.null(T)) {
      dim_outcome_array <- dim(old_outcome_array)
      name_outcome_array <- dimnames(old_outcome_array)
      levels <- (1:dim_outcome_array[3])[T]
      if (type == "all")
        contrast_levels <- (function(x) x[x[,1] > x[,2], ])(expand.grid(levels, levels))
      if (type == "sequential" || (type == "cycle" && length(levels) == 2))
        contrast_levels <- data.frame(Var1 = levels[-1], Var2 = levels[-length(levels)])
      if (type == "cycle" && length(levels) > 2) {
        contrast_levels <- data.frame(Var1 = c(levels[-1], levels[1]), Var2 = levels)
      }
      dim_outcome_array[3] <- nrow(contrast_levels)
      name_outcome_array[[3]] <- apply(contrast_levels, 1, function(x) {
        paste0("{", name_outcome_array[[3]][x[1]], "}-{", name_outcome_array[[3]][x[2]], "}")
      })
      new_outcome_array <- array(NA, dim = dim_outcome_array, dimnames = name_outcome_array)
      for (i in 1:nrow(contrast_levels)) {
        new_outcome_array[, , i, ] <- old_outcome_array[, , contrast_levels[i, 1] ,] - 
          old_outcome_array[, , contrast_levels[i, 2], ]
      }
      old_outcome_array <- new_outcome_array
    }
  }
  return (structure(
    list(
      outcome_array = old_outcome_array,
      is.survival = outcome$is.survival,
      time_points = outcome$time_points
    ),
    class = c("PSContrast", "PSOutcome")
  ))
}

#' @export 
print.PSContrast <- function(x, ...) {
  if (!x$is.survival) {
    cat("Non-survival PSContrast Object (", dim(x$outcome_array)[1],
        " strata, ", dim(x$outcome_array)[2], 
        " treatment arms, ", dim(x$outcome_array)[3], " iterations)\n", sep = '')
  } else {
    cat("Survival PSContrast Object (", dim(x$outcome_array)[1],
        " strata, ", dim(x$outcome_array)[2], 
        " treatment arms, ", dim(x$outcome_array)[4], " iterations)\n", sep = '')
    cat("Evaluated at ", dim(x$outcome_array)[3], " time points.\n", sep = '')
  }
}
  