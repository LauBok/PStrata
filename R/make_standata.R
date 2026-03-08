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
make_standata <- function(PSobject) {
  stopifnot(
    "Rows of the data frames provided for S.formula and Y.formula must match" =
      nrow(PSobject$S.formula$data) == nrow(PSobject$Y.formula$data)
  )

  D_names <- PSobject$S.formula$response_names[-1]
  D_data <- dplyr::select(PSobject$S.formula$data, dplyr::all_of(D_names))
  D_int <- apply(D_data, 1, function(x) {
    t(c(x, 0)) %*% c(1, cumprod(PSobject$strata_info$max_postrand_level + 1))
  })

  Z <- dplyr::pull(PSobject$S.formula$data, PSobject$S.formula$response_names[1])
  Z_int <- encode_treatment(Z, PSobject$strata_info$num_treatment)

  df <- list(
    N  = nrow(PSobject$S.formula$data),
    PS = ncol(PSobject$S.formula$fixed_eff_matrix),
    PG = ncol(PSobject$Y.formula$fixed_eff_matrix),
    Z  = Z_int,
    D  = D_int,
    Y  = dplyr::pull(PSobject$Y.formula$data, PSobject$Y.formula$response_names[1]),
    XS = PSobject$S.formula$fixed_eff_matrix,
    XG = PSobject$Y.formula$fixed_eff_matrix
  )

  # Survival-specific data
  if (PSobject$is.survival) {
    df$delta <- dplyr::pull(PSobject$Y.formula$data, PSobject$Y.formula$response_names[2])
    df$time <- if (length(PSobject$survival.time.points) == 1) {
      seq(0, max(df$Y) * 0.9, length.out = PSobject$survival.time.points)
    } else {
      PSobject$survival.time.points
    }
    df$T <- length(df$time)
  }

  # Random effects for S-model
  df <- append_re_data(df, PSobject$S.formula$random_eff_list, "S")
  # Random effects for Y-model
  df <- append_re_data(df, PSobject$Y.formula$random_eff_list, "G")

  df
}

#' Encode treatment variable as integer (0-indexed)
#' @noRd
encode_treatment <- function(Z, num_treatment) {
  if (is.factor(Z)) {
    as.integer(Z) - 1
  } else if (all(sapply(Z, function(x) x %in% 0:(num_treatment - 1)))) {
    Z
  } else {
    as.integer(as.factor(Z)) - 1
  }
}

#' Append random effect data to the Stan data list
#' @param df The Stan data list
#' @param re_list Random effect list from PSFormula
#' @param prefix "S" or "G"
#' @noRd
append_re_data <- function(df, re_list, prefix) {
  if (is.null(re_list)) return(df)
  for (i in seq_along(re_list)) {
    df[[paste0("P", prefix, "_RE_", i)]] <- length(re_list[[i]]$terms)
    df[[paste0("N", prefix, "_RE_", i)]] <- length(attr(re_list[[i]]$factors, "levels"))
    df[[paste0("X", prefix, "_RE_", i)]] <- re_list[[i]]$matrix
  }
  df
}
