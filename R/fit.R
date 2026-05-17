# ===========================================================================
# fit: Estimate a PStrataModel with data via Stan MCMC
# ===========================================================================

#' Fit a Principal Stratification Model
#'
#' Estimate a principal stratification model by running MCMC via Stan.
#' The model specification (a \code{\link{PStrataModel}}) is combined with
#' observed data to produce posterior samples.
#'
#' @param model A \code{\link{PStrataModel}} object.
#' @param data A data frame containing all variables referenced in the model.
#' @param chains,iter,warmup,cores,seed MCMC settings passed to
#'   \code{\link[rstan]{stan}}.
#' @param .debug Logical. If TRUE, read Stan helper files from \code{inst/}
#'   rather than installed package data. Useful during development.
#' @param ... Additional arguments passed to \code{\link[rstan]{stan}}.
#'
#' @return An object of class \code{PStrataFit}
#'   (or \code{PStrataFitSurvival}).
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
#' summary(ps_fit)
#' diagnostics(ps_fit)
#' plot(ps_fit)
#' cat(ps_fit$stancode)
#' }
#' @export
fit <- function(model, ...) UseMethod("fit")

#' @rdname fit
#' @export
fit.PStrataModel <- function(model, data, chains = 4, iter = 2000,
                              warmup = floor(iter / 2), cores = 1,
                              seed = NULL, .debug = FALSE, ...) {
  cl <- match.call()
  stopifnot(is.data.frame(data))

  # 1. Build design matrices and extract columns from data
  dm <- build_design_matrices(model, data)

  # 2. Encode strata definition as integer matrix
  strata_enc <- encode_strata(model$strata_def, model$num_treatment)

  # 3. Build SZDG mapping table
  SZDG <- build_SZDG(
    strata_enc$strata_matrix, model$exclusion_restriction,
    model$num_strata, model$num_treatment
  )

  # 4. Assemble Stan data list
  standata <- assemble_standata(model, dm, strata_enc, SZDG)

  # 5. Generate Stan code via bridge to existing make_stancode
  bridge <- make_stancode_bridge(model, dm, SZDG)
  stancode <- make_stancode(bridge, debug = .debug)

  # 6. Compile and run MCMC
  stanmodel <- rstan::stan_model(model_code = stancode)
  sampling_args <- list(
    object = stanmodel, data = standata,
    chains = chains, iter = iter, warmup = warmup, cores = cores
  )
  if (!is.null(seed)) sampling_args$seed <- seed
  stanfit <- do.call(rstan::sampling, c(sampling_args, list(...)))

  # 7. Return PStrataFit
  res <- list(
    call               = cl,
    model              = model,
    data               = data,
    stanfit            = stanfit,
    stancode           = stancode,
    standata           = standata,
    SZDG_table         = SZDG,
    strata_matrix      = strata_enc$strata_matrix,
    max_postrand_level = strata_enc$max_postrand_level,
    design_matrices    = list(XS = dm$XS, XG = dm$XG),
    time_points        = if (model$is_survival) standata$time else NULL
  )

  class(res) <- if (model$is_survival) {
    c("PStrataFitSurvival", "PStrataFit")
  } else {
    "PStrataFit"
  }

  res
}

# ===========================================================================
# Design matrix construction
# ===========================================================================

#' Build design matrices and extract response columns from data
#' @noRd
build_design_matrices <- function(model, data) {
  # Validate LHS columns exist
  needed <- c(model$treatment_name, model$postrand_names, model$outcome_name)
  if (!is.null(model$event_name)) needed <- c(needed, model$event_name)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0)
    stop("Columns not found in data: ", paste(missing_cols, collapse = ", "))

  # Validate RE grouping variables exist
  re_groups <- c(
    vapply(model$s_random_effects %||% list(), `[[`, character(1), "group"),
    vapply(model$y_random_effects %||% list(), `[[`, character(1), "group")
  )
  missing_re <- setdiff(re_groups, names(data))
  if (length(missing_re) > 0)
    stop("RE grouping variables not found in data: ",
         paste(missing_re, collapse = ", "))

  # Fixed-effect design matrices (strip response and random effects)
  s_terms <- stats::delete.response(
    stats::terms(reformulas::nobars(model$stratum_formula))
  )
  y_terms <- stats::delete.response(
    stats::terms(reformulas::nobars(model$outcome_formula))
  )
  XS <- stats::model.matrix(s_terms, data)
  XG <- stats::model.matrix(y_terms, data)

  # Extract response columns
  Z      <- data[[model$treatment_name]]
  D_cols <- lapply(model$postrand_names, function(d) data[[d]])
  Y      <- data[[model$outcome_name]]
  delta  <- if (!is.null(model$event_name)) data[[model$event_name]] else NULL

  # Random effect design matrices (via lme4)
  s_re <- parse_re_from_data(model$stratum_formula, model$s_random_effects, data)
  y_re <- parse_re_from_data(model$outcome_formula, model$y_random_effects, data)

  list(
    XS = XS, XG = XG,
    Z = Z, D_cols = D_cols, Y = Y, delta = delta,
    s_re = s_re, y_re = y_re
  )
}

#' Parse random effect design matrices from actual data using lme4
#' @noRd
parse_re_from_data <- function(formula, re_info, data) {
  if (is.null(re_info) || length(re_info) == 0) return(NULL)

  # Construct dummy-response formula for lme4
  rhs_str <- paste(deparse(formula[[3]], width.cutoff = 500), collapse = "")
  data$.y <- 0
  dummy_f <- stats::as.formula(paste(".y ~", rhs_str))

  lf <- lme4::lFormula(dummy_f, data)

  mapply(
    function(mat, cnms, fac) {
      list(matrix = t(as.matrix(mat)), terms = cnms, factors = fac)
    },
    lf$reTrms$Ztlist,
    lf$reTrms$cnms,
    lf$reTrms$flist,
    SIMPLIFY = FALSE
  )
}

# ===========================================================================
# Strata encoding
# ===========================================================================

#' Encode strata_def into integer strata_matrix for Stan
#'
#' Converts the list-of-lists strata definition into an integer matrix
#' where rows = strata, cols = treatment arms, values = encoded D integer.
#' Multiple D variables are encoded via polynomial: d1 + d2*(m1+1) + ...
#'
#' @noRd
encode_strata <- function(strata_def, num_treatment) {
  num_postrand <- length(strata_def[[1]])

  # Compute max level per postrand variable
  max_level <- integer(num_postrand)
  for (s in strata_def) {
    for (j in seq_len(num_postrand)) {
      max_level[j] <- max(max_level[j], max(s[[j]]))
    }
  }

  # Polynomial encoding bases
  bases <- c(1L, cumprod(max_level + 1L))

  mat <- matrix(0L, nrow = length(strata_def), ncol = num_treatment)
  for (i in seq_along(strata_def)) {
    for (z in seq_len(num_treatment)) {
      d_vals <- vapply(strata_def[[i]], function(v) v[z], integer(1))
      mat[i, z] <- as.integer(sum(c(d_vals, 0L) * bases))
    }
  }

  list(strata_matrix = mat, max_postrand_level = max_level)
}

#' Encode observed D columns from data into a single integer vector
#' @noRd
encode_D_data <- function(D_cols, max_postrand_level) {
  if (length(D_cols) == 1) return(as.integer(D_cols[[1]]))
  D_mat <- do.call(cbind, D_cols)
  bases <- c(1L, cumprod(max_postrand_level + 1L))
  apply(D_mat, 1, function(x) as.integer(sum(c(x, 0L) * bases)))
}

#' Encode treatment values to 0-based integers
#' @noRd
encode_Z_data <- function(Z, num_treatment) {
  if (is.factor(Z)) return(as.integer(Z) - 1L)
  z_int <- as.integer(Z)
  if (all(z_int %in% 0L:(num_treatment - 1L))) return(z_int)
  warning("Treatment variable auto-converted to 0-based integers.")
  as.integer(as.factor(Z)) - 1L
}

# ===========================================================================
# SZDG table construction
# ===========================================================================

#' Build the (Stratum, Treatment, D-encoded, outcome Group) mapping table
#'
#' Each row maps a combination of stratum and treatment arm to the observed
#' D value and the outcome parameter group. Under exclusion restriction (ER),
#' a stratum with identical D values across treatment arms shares one group.
#'
#' @noRd
build_SZDG <- function(strata_matrix, er_list, num_strata, num_treatment) {
  SZDG <- matrix(nrow = 0L, ncol = 4L)
  colnames(SZDG) <- c("S", "Z", "D", "G")
  G_id <- 0L

  for (s in seq_len(num_strata)) {
    S_id <- s - 1L
    for (z in seq_len(num_treatment)) {
      Z_id <- z - 1L
      D_id <- strata_matrix[s, z]

      # Under ER, reuse group if same stratum & D value already seen
      tmp_G <- NULL
      if (er_list[s] && nrow(SZDG) > 0) {
        match_idx <- which(SZDG[, "S"] == S_id & SZDG[, "D"] == D_id)
        if (length(match_idx) > 0) tmp_G <- SZDG[match_idx[1], "G"]
      }

      if (is.null(tmp_G)) {
        G_id <- G_id + 1L
        tmp_G <- G_id
      }
      SZDG <- rbind(SZDG, c(S_id, Z_id, D_id, tmp_G))
    }
  }

  SZDG
}

# ===========================================================================
# Stan data assembly
# ===========================================================================

#' Assemble the named list of data for Stan
#' @noRd
assemble_standata <- function(model, dm, strata_enc, SZDG) {
  Z_int <- encode_Z_data(dm$Z, model$num_treatment)
  D_int <- encode_D_data(dm$D_cols, strata_enc$max_postrand_level)

  sd <- list(
    N  = nrow(dm$XS),
    PS = ncol(dm$XS),
    PG = ncol(dm$XG),
    Z  = Z_int,
    D  = D_int,
    Y  = dm$Y,
    XS = dm$XS,
    XG = dm$XG
  )

  # Survival-specific fields
  if (model$is_survival) {
    sd$delta <- as.integer(dm$delta)
    stp <- model$survival_time_points
    sd$time <- if (length(stp) == 1) {
      seq(0, max(sd$Y) * 0.9, length.out = stp)
    } else {
      stp
    }
    sd[["T"]] <- length(sd$time)
  }

  # Random effect matrices
  sd <- append_re_standata(sd, dm$s_re, "S")
  sd <- append_re_standata(sd, dm$y_re, "G")

  sd
}

#' Append random effect data to the Stan data list
#' @noRd
append_re_standata <- function(sd, re_list, prefix) {
  if (is.null(re_list)) return(sd)
  for (i in seq_along(re_list)) {
    sd[[paste0("P", prefix, "_RE_", i)]] <- length(re_list[[i]]$terms)
    sd[[paste0("N", prefix, "_RE_", i)]] <- length(levels(re_list[[i]]$factors))
    sd[[paste0("X", prefix, "_RE_", i)]] <- re_list[[i]]$matrix
  }
  sd
}

# ===========================================================================
# Bridge to existing make_stancode
# ===========================================================================

#' Create a minimal object that the existing make_stancode / build_stan_context
#' can consume, mapping PStrataModel fields to the PSObject-like interface.
#' @noRd
make_stancode_bridge <- function(model, dm, SZDG) {
  safe <- function(x) if (is.null(x)) prior_flat() else x

  list(
    SZDG_table        = SZDG,
    Y.family          = model$outcome_family,
    S.formula         = list(random_eff_list = dm$s_re),
    Y.formula         = list(random_eff_list = dm$y_re),
    prior_intercept   = model$priors$intercept,
    prior_coefficient = model$priors$coefficient,
    prior_sigma       = safe(model$priors$sigma),
    prior_alpha       = safe(model$priors$alpha),
    prior_lambda      = safe(model$priors$lambda),
    prior_theta       = safe(model$priors$theta)
  )
}

# ===========================================================================
# Null-default operator (if not already available)
# ===========================================================================

if (!exists("%||%")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# ===========================================================================
# PStrataFit: print
# ===========================================================================

#' @export
print.PStrataFit <- function(x, ...) {
  m <- x$model
  samples <- rstan::extract(x$stanfit)

  cat("PStrataFit\n")
  cat("  Family:", m$outcome_family$family,
      "(", m$outcome_family$link, ")\n")

  # Strata with ER markers
  strata_strs <- vapply(seq_along(m$strata_names), function(i) {
    if (m$exclusion_restriction[i]) paste0(m$strata_names[i], " (ER)")
    else m$strata_names[i]
  }, character(1))
  cat("  Strata:", paste(strata_strs, collapse = ", "), "\n")
  cat("  Outcome groups:", paste(m$outcome_groups$labels, collapse = ", "), "\n")

  cat("  N:", x$standata$N, "observations\n")

  si <- x$stanfit@sim
  cat("  MCMC:", si$chains, "chains,", si$iter, "iter",
      "(", si$warmup, "warmup )\n")

  # Stratum proportions
  strata_prob <- colMeans(samples$strata_prob)
  names(strata_prob) <- m$strata_names
  cat("\nEstimated stratum proportions:\n")
  print(round(strata_prob, 4))

  # Mean outcome per group (non-survival only)
  if ("mean_effect" %in% names(samples)) {
    mean_eff <- colMeans(as.matrix(samples$mean_effect))
    names(mean_eff) <- m$outcome_groups$labels
    cat("\nEstimated mean outcome per group:\n")
    print(round(mean_eff, 4))
  }

  cat("\nUse summary() for posterior intervals.\n")
  cat("Use diagnostics() to check MCMC convergence.\n")
  invisible(x)
}

#' @export
print.PStrataFitSurvival <- function(x, ...) {
  NextMethod()
  tp <- x$time_points
  if (!is.null(tp) && length(tp) > 1) {
    cat("Survival evaluated at", length(tp), "time points in [",
        round(min(tp), 2), ",", round(max(tp), 2), "].\n")
  }
}

# ===========================================================================
# PStrataFit: summary
# ===========================================================================

#' @export
summary.PStrataFit <- function(object,
                                probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                ...) {
  samples <- rstan::extract(object$stanfit)
  m <- object$model
  prob_names <- paste0(probs * 100, "%")

  # Stratum proportions
  sp <- samples$strata_prob
  sp_summary <- t(apply(sp, 2, function(x) {
    c(mean(x), stats::sd(x), stats::quantile(x, probs))
  }))
  rownames(sp_summary) <- m$strata_names
  colnames(sp_summary) <- c("mean", "sd", prob_names)

  # Mean outcome per group (non-survival)
  mean_effect_summary <- NULL
  if ("mean_effect" %in% names(samples)) {
    me <- as.matrix(samples$mean_effect)
    mean_effect_summary <- t(apply(me, 2, function(x) {
      c(mean(x), stats::sd(x), stats::quantile(x, probs))
    }))
    rownames(mean_effect_summary) <- m$outcome_groups$labels
    colnames(mean_effect_summary) <- c("mean", "sd", prob_names)
  }

  # beta_S (reference = first stratum)
  beta_S_summary <- summarize_coef_array(
    samples$beta_S,
    group_names = m$strata_names[-1],
    coef_names  = colnames(object$design_matrices$XS),
    probs       = probs
  )

  # beta_G
  beta_G_summary <- summarize_coef_array(
    samples$beta_G,
    group_names = m$outcome_groups$labels,
    coef_names  = colnames(object$design_matrices$XG),
    probs       = probs
  )

  # Extra parameters (sigma, theta, etc.)
  extra_summary <- list()
  for (ep in m$family_info$extra_params) {
    if (ep$name %in% names(samples)) {
      vals <- as.matrix(samples[[ep$name]])
      ep_sum <- t(apply(vals, 2, function(x) {
        c(mean(x), stats::sd(x), stats::quantile(x, probs))
      }))
      rownames(ep_sum) <- m$outcome_groups$labels
      colnames(ep_sum) <- c("mean", "sd", prob_names)
      extra_summary[[ep$name]] <- ep_sum
    }
  }

  res <- list(
    call          = object$call,
    model         = m,
    strata_prob   = sp_summary,
    mean_effect   = mean_effect_summary,
    beta_S        = beta_S_summary,
    beta_G        = beta_G_summary,
    extra_params  = extra_summary
  )
  class(res) <- "summary.PStrataFit"
  res
}

#' Summarize a 3D posterior coefficient array (iterations x groups x covariates)
#' @noRd
summarize_coef_array <- function(arr, group_names, coef_names, probs) {
  n_groups <- dim(arr)[2]
  n_coefs  <- dim(arr)[3]
  prob_names <- paste0(probs * 100, "%")

  rows   <- vector("list", n_groups * n_coefs)
  rnames <- character(n_groups * n_coefs)
  k <- 0L
  for (g in seq_len(n_groups)) {
    for (p in seq_len(n_coefs)) {
      k <- k + 1L
      x <- arr[, g, p]
      rows[[k]] <- c(mean(x), stats::sd(x), stats::quantile(x, probs))
      rnames[k] <- paste0(group_names[g], ":", coef_names[p])
    }
  }

  mat <- do.call(rbind, rows)
  rownames(mat) <- rnames
  colnames(mat) <- c("mean", "sd", prob_names)
  mat
}

#' @export
print.summary.PStrataFit <- function(x, digits = 4, ...) {
  cat("PStrataFit Summary\n")
  cat("====================\n\n")

  cat("Call:", deparse(x$call, width.cutoff = 200), "\n\n")

  cat("Stratum proportions:\n")
  print(round(x$strata_prob, digits))

  if (!is.null(x$mean_effect)) {
    cat("\nMean outcome per group:\n")
    print(round(x$mean_effect, digits))
  }

  cat("\nStratum model coefficients (beta_S, reference:",
      x$model$strata_names[1], "):\n")
  print(round(x$beta_S, digits))

  cat("\nOutcome model coefficients (beta_G):\n")
  print(round(x$beta_G, digits))

  for (nm in names(x$extra_params)) {
    cat("\n", nm, ":\n", sep = "")
    print(round(x$extra_params[[nm]], digits))
  }

  invisible(x)
}

# ===========================================================================
# PStrataFit: diagnostics
# ===========================================================================

#' MCMC Convergence Diagnostics
#'
#' Reports the rank-normalized split-Rhat together with the bulk and tail
#' effective sample sizes (Vehtari et al., 2021), as computed by
#' \code{\link[rstan]{monitor}}.
#'
#' @param object A \code{PStrataFit} object.
#' @param pars Character vector of parameter names. If NULL, key parameters
#'   are selected automatically.
#' @param ... Currently unused.
#' @return Invisibly returns the \code{rstan::monitor} summary.
#' @export
diagnostics <- function(object, ...) UseMethod("diagnostics")

#' @rdname diagnostics
#' @export
diagnostics.PStrataFit <- function(object, pars = NULL, ...) {
  if (is.null(pars)) {
    pars <- c("beta_S", "beta_G", "strata_prob")
    all_pars <- names(rstan::extract(object$stanfit))
    extra <- intersect(c("sigma", "alpha", "lambda", "theta"), all_pars)
    pars <- c(pars, extra)
  }

  sims <- rstan::extract(object$stanfit, pars = pars, permuted = FALSE)
  mon  <- as.data.frame(
    suppressWarnings(rstan::monitor(sims, warmup = 0, print = FALSE))
  )
  rhat <- stats::setNames(mon[, "Rhat"], rownames(mon))
  bulk <- mon[, "Bulk_ESS"]
  tail <- mon[, "Tail_ESS"]

  cat("MCMC Diagnostics\n")
  cat("----------------\n")
  cat("Rhat (rank-normalized split-Rhat) - max:",
      round(max(rhat, na.rm = TRUE), 4),
      " min:", round(min(rhat, na.rm = TRUE), 4), "\n")
  cat("Bulk-ESS - min:", round(min(bulk, na.rm = TRUE)), "\n")
  cat("Tail-ESS - min:", round(min(tail, na.rm = TRUE)), "\n")

  if (any(rhat > 1.01, na.rm = TRUE)) {
    bad <- names(rhat[rhat > 1.01])
    cat("WARNING: Rhat > 1.01 for:",
        paste(utils::head(bad, 10), collapse = ", "),
        if (length(bad) > 10) "...", "\n")
  } else {
    cat("All Rhat values <= 1.01 (good).\n")
  }

  invisible(mon)
}

# ===========================================================================
# PStrataFit: accessors
# ===========================================================================

#' @export
coef.PStrataFit <- function(object, ...) {
  samples <- rstan::extract(object$stanfit)
  result <- list()

  if ("beta_S" %in% names(samples)) {
    beta_S <- apply(samples$beta_S, c(2, 3), mean)
    rownames(beta_S) <- object$model$strata_names[-1]
    colnames(beta_S) <- colnames(object$design_matrices$XS)
    result$beta_S <- beta_S
  }

  if ("beta_G" %in% names(samples)) {
    beta_G <- apply(samples$beta_G, c(2, 3), mean)
    rownames(beta_G) <- object$model$outcome_groups$labels
    colnames(beta_G) <- colnames(object$design_matrices$XG)
    result$beta_G <- beta_G
  }

  for (p in c("sigma", "alpha", "lambda", "theta")) {
    if (p %in% names(samples)) {
      vals <- colMeans(as.matrix(samples[[p]]))
      names(vals) <- object$model$outcome_groups$labels
      result[[p]] <- vals
    }
  }

  result
}

#' @export
nobs.PStrataFit <- function(object, ...) object$standata$N

#' @export
formula.PStrataFit <- function(x, ...) {
  list(S = x$model$stratum_formula, Y = x$model$outcome_formula)
}

#' Extract the Generated Stan Code
#'
#' @param object A \code{PStrataFit} object.
#' @param ... Currently unused.
#' @return The Stan code as a character string.
#' @export
stancode.PStrataFit <- function(object, ...) object$stancode

# ===========================================================================
# PStrataFit: plot
# ===========================================================================

#' Plot PStrataFit Objects
#'
#' By default, shows violin plots of posterior strata proportions
#' and potential outcomes (non-survival).
#' Use \code{what = "trace"} or \code{what = "dens"} for MCMC diagnostics.
#'
#' @param x A \code{PStrataFit} object.
#' @param what \code{"default"} for strata proportions + potential outcomes,
#'   \code{"trace"} for traceplots, \code{"dens"} for density plots.
#' @param pars Parameter names for trace/dens plots (ignored for default).
#' @param ... Additional arguments passed to plotting functions.
#' @return A ggplot object.
#' @export
plot.PStrataFit <- function(x, what = c("default", "trace", "dens"),
                             pars = "strata_prob", ...) {
  what <- match.arg(what)
  if (what == "trace") return(rstan::traceplot(x$stanfit, pars = pars, ...))
  if (what == "dens") return(rstan::stan_dens(x$stanfit, pars = pars, ...))

  m <- x$model
  samples <- rstan::extract(x$stanfit)

  # --- Strata proportions: [iter, S] → long ---
  sp <- samples$strata_prob
  sp_df <- data.frame(
    label = paste0("P(", rep(m$strata_names, each = nrow(sp)), ")"),
    value = as.vector(sp),
    stringsAsFactors = FALSE
  )

  # --- Potential outcomes (non-survival) ---
  if (!m$is_survival && "mean_effect" %in% names(samples)) {
    est <- estimate(x)
    arr <- est$outcome_array  # [S, Z, iter]
    dims <- dim(arr)
    dn <- dimnames(arr)

    po_grid <- expand.grid(
      iter = seq_len(dims[3]),
      Z    = dn[[2]],
      S    = dn[[1]],
      stringsAsFactors = FALSE
    )
    po_df <- data.frame(
      label = paste0("E[Y|", po_grid$S, ",", po_grid$Z, "]"),
      value = as.vector(aperm(arr, c(3, 2, 1))),
      stringsAsFactors = FALSE
    )
    plot_df <- rbind(sp_df, po_df)
  } else {
    plot_df <- sp_df
  }
  plot_df$label <- factor(plot_df$label, levels = unique(plot_df$label))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = value)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.5) +
    ggplot2::facet_wrap(~ label, scales = "free") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y  = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
}
