# Internal utility functions for PStrata
# These are not exported and are used by multiple files.

# --- Array helpers ---

#' Subset an array along a specific dimension
#' @param arr An array
#' @param dim_idx Integer, the dimension to subset along
#' @param idx Index (numeric, logical, or character) to select
#' @return A subsetted array with drop = FALSE
#' @noRd
slice_array <- function(arr, dim_idx, idx) {
  indices <- rep(list(TRUE), length(dim(arr)))
  indices[[dim_idx]] <- idx
  do.call(`[`, c(list(arr), indices, list(drop = FALSE)))
}

#' Assign values to a slice of an array along a specific dimension
#' @noRd
assign_slice <- function(arr, dim_idx, idx, value) {
  indices <- rep(list(TRUE), length(dim(arr)))
  indices[[dim_idx]] <- idx
  do.call(`[<-`, c(list(arr), indices, list(value = value)))
}

# --- Contrast computation ---

#' Compute pairs of indices for contrast computation
#' @param levels Integer vector of dimension indices to contrast
#' @param type One of "all", "sequential", "cycle"
#' @return A data.frame with columns Var1, Var2 (each row is a pair to subtract)
#' @noRd
compute_contrast_pairs <- function(levels, type = c("all", "sequential", "cycle")) {
  type <- match.arg(type)
  if (type == "all") {
    pairs <- expand.grid(Var1 = levels, Var2 = levels)
    pairs[pairs[, 1] > pairs[, 2], , drop = FALSE]
  } else if (type == "sequential" || (type == "cycle" && length(levels) == 2)) {
    data.frame(Var1 = levels[-1], Var2 = levels[-length(levels)])
  } else { # cycle with > 2 levels
    data.frame(Var1 = c(levels[-1], levels[1]), Var2 = levels)
  }
}

#' Compute contrasts along a single dimension of an array
#'
#' Replaces the repetitive contrast logic in PSContrast for each dimension
#' and for survival/non-survival cases.
#'
#' @param arr An array (3D for non-survival, 4D for survival)
#' @param dim_idx Which dimension to take contrasts along (1=S, 2=Z, 3=T)
#' @param selection TRUE (all levels) or a logical/numeric index vector
#' @param type One of "all", "sequential", "cycle"
#' @return A new array with contrasts computed along the specified dimension
#' @noRd
contrast_along_dim <- function(arr, dim_idx, selection, type) {
  dim_arr <- dim(arr)
  name_arr <- dimnames(arr)

  if (is.logical(selection) && length(selection) == 1 && selection) {
    levels <- seq_len(dim_arr[dim_idx])
  } else {
    levels <- (seq_len(dim_arr[dim_idx]))[selection]
  }

  pairs <- compute_contrast_pairs(levels, type)

  new_dim <- dim_arr
  new_dim[dim_idx] <- nrow(pairs)

  new_names <- name_arr
  new_names[[dim_idx]] <- apply(pairs, 1, function(x) {
    paste0("{", name_arr[[dim_idx]][x[1]], "}-{", name_arr[[dim_idx]][x[2]], "}")
  })

  new_arr <- array(NA, dim = new_dim, dimnames = new_names)

  for (i in seq_len(nrow(pairs))) {
    diff <- slice_array(arr, dim_idx, pairs[i, 1]) -
            slice_array(arr, dim_idx, pairs[i, 2])
    new_arr <- assign_slice(new_arr, dim_idx, i, diff)
  }

  new_arr
}

# --- Summary helpers ---

#' Compute summary statistics (mean, sd, quantiles) over the last dimension
#'
#' The last dimension of the array is treated as the iteration/sample dimension.
#' Returns an array where the last dimension is the 7 summary statistics.
#'
#' @param arr An array where the last dimension is iterations
#' @return An array with the last dimension replaced by 7 summary stats
#' @noRd
summarize_last_dim <- function(arr) {
  ndim <- length(dim(arr))
  margin <- rev(seq_len(ndim - 1))
  summary_raw <- apply(arr, margin, function(x) {
    c(mean(x), stats::sd(x),
      stats::quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))
  })
  aperm(summary_raw, rev(seq_len(length(dim(summary_raw)))))
}

#' Standard summary statistic names
#' @noRd
summary_stat_names <- function() {
  c("mean", "sd", "2.5%", "25%", "median", "75%", "97.5%")
}

#' Convert a summary array to a matrix with labeled rows
#'
#' Flattens all dimensions except the last (statistics) into rows.
#' Row labels are colon-separated combinations of dimension names.
#' The ordering matches the original: first dim varies slowest, last data dim fastest.
#'
#' @param summary_array Array where last dimension is summary statistics
#' @return A matrix with rows labeled "dim1:dim2:..." and columns as stat names
#' @noRd
summary_array_to_matrix <- function(summary_array) {
  dims <- dim(summary_array)
  ndim <- length(dims)
  n_stats <- dims[ndim]
  data_dims <- dims[-ndim]
  n_rows <- prod(data_dims)

  mat <- matrix(NA, nrow = n_rows, ncol = n_stats)

  # Build index grid: last data dim varies fastest (innermost loop)
  idx_grid <- do.call(expand.grid, rev(lapply(data_dims, seq_len)))
  idx_grid <- idx_grid[, rev(seq_len(ncol(idx_grid))), drop = FALSE]

  for (i in seq_len(n_rows)) {
    idx <- c(as.list(idx_grid[i, ]), list(TRUE))
    mat[i, ] <- do.call(`[`, c(list(summary_array), idx, list(drop = TRUE)))
  }

  # Row names from dimension name combinations
  dim_names <- dimnames(summary_array)
  name_grid <- do.call(expand.grid, rev(lapply(dim_names[-ndim], identity)))
  name_grid <- name_grid[, rev(seq_len(ncol(name_grid))), drop = FALSE]
  rownames(mat) <- apply(name_grid, 1, paste, collapse = ":")
  colnames(mat) <- dim_names[[ndim]]

  mat
}

#' Convert a summary matrix to a data.frame with dimension columns
#'
#' Splits row names by ":" and prepends them as named columns.
#'
#' @param mat A summary matrix from summary_array_to_matrix
#' @param col_names Character vector of column names for the dimension columns
#' @param time_points Optional numeric vector; if provided, integer T values are
#'   converted to the corresponding time point values
#' @return A data.frame with dimension columns prepended
#' @noRd
summary_matrix_to_df <- function(mat, col_names, time_points = NULL) {
  df <- as.data.frame(mat)
  parts <- purrr::transpose(stringr::str_split(rownames(df), ":"))

  prefix_df <- as.data.frame(
    stats::setNames(
      lapply(seq_along(col_names), function(i) unlist(parts[[i]])),
      col_names
    ),
    stringsAsFactors = FALSE
  )

  # Convert T column to numeric time points if applicable
  if ("T" %in% col_names && !is.null(time_points)) {
    t_col <- prefix_df$T
    t_int <- suppressWarnings(as.integer(t_col))
    if (!any(is.na(t_int))) {
      prefix_df$T <- time_points[t_int]
    }
  }

  cbind(prefix_df, df)
}
