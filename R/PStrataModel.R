# ===========================================================================
# PStrataModel: Model specification for principal stratification
# ===========================================================================

# --- Family registry ---
# Maps family + link to outcome type, Stan functions, and display equations

.family_registry <- list(
  gaussian = list(
    outcome_type = "real",
    stan_func = "normal_lpdf",
    extra_params = list(list(name = "sigma", type = "positive")),
    links = list(
      identity = list(
        stan_link = "identity_fn",
        eq = "Y | g, X ~ Normal(eta_g, sigma_g)"
      ),
      log = list(
        stan_link = "exp",
        eq = "Y | g, X ~ Normal(exp(eta_g), sigma_g)"
      ),
      inverse = list(
        stan_link = "inv",
        eq = "Y | g, X ~ Normal(1 / eta_g, sigma_g)"
      )
    )
  ),
  binomial = list(
    outcome_type = "binary",
    stan_func = "bernoulli_logit_lpmf",
    extra_params = list(),
    links = list(
      logit = list(
        stan_link = "identity_fn",
        eq = "P(Y=1 | g, X) = logistic(eta_g)"
      ),
      probit = list(
        stan_link = "inv_Phi",
        eq = "P(Y=1 | g, X) = Phi(eta_g)"
      ),
      cauchit = list(
        stan_link = "cauchy_cdf",
        eq = "P(Y=1 | g, X) = cauchy_cdf(eta_g)"
      ),
      log = list(
        stan_link = "exp",
        eq = "P(Y=1 | g, X) = exp(eta_g)"
      ),
      cloglog = list(
        stan_link = "inv_cloglog",
        eq = "P(Y=1 | g, X) = 1 - exp(-exp(eta_g))"
      )
    )
  ),
  Gamma = list(
    outcome_type = "positive",
    stan_func = "gamma_lpdf",
    extra_params = list(list(name = "alpha", type = "positive")),
    links = list(
      inverse = list(
        stan_link = "inv",
        eq = "Y | g, X ~ Gamma(alpha_g, alpha_g * eta_g)"
      ),
      identity = list(
        stan_link = "identity_fn",
        eq = "Y | g, X ~ Gamma(alpha_g, alpha_g / eta_g)"
      ),
      log = list(
        stan_link = "exp",
        eq = "Y | g, X ~ Gamma(alpha_g, alpha_g / exp(eta_g))"
      )
    )
  ),
  poisson = list(
    outcome_type = "count",
    stan_func = "poisson_lpmf",
    extra_params = list(),
    links = list(
      log = list(
        stan_link = "exp",
        eq = "Y | g, X ~ Poisson(exp(eta_g))"
      ),
      identity = list(
        stan_link = "identity_fn",
        eq = "Y | g, X ~ Poisson(eta_g)"
      ),
      sqrt = list(
        stan_link = "square",
        eq = "Y | g, X ~ Poisson(eta_g^2)"
      )
    )
  ),
  inverse.gaussian = list(
    outcome_type = "positive",
    stan_func = "inv_gaussian_lpdf",
    extra_params = list(list(name = "lambda", type = "positive")),
    links = list(
      "1/mu^2" = list(
        stan_link = "inv_sqrt",
        eq = "Y | g, X ~ InvGaussian(1 / sqrt(eta_g), lambda_g)"
      ),
      inverse = list(
        stan_link = "inv",
        eq = "Y | g, X ~ InvGaussian(1 / eta_g, lambda_g)"
      ),
      identity = list(
        stan_link = "identity_fn",
        eq = "Y | g, X ~ InvGaussian(eta_g, lambda_g)"
      ),
      log = list(
        stan_link = "exp",
        eq = "Y | g, X ~ InvGaussian(exp(eta_g), lambda_g)"
      )
    )
  ),
  survival_Cox = list(
    outcome_type = "survival",
    stan_func = "weibull_cox_lpdf",
    extra_params = list(list(name = "theta", type = "positive")),
    links = list(
      identity = list(
        stan_link = "identity_fn",
        eq = "h(t | g, X) = theta_g * t^(theta_g - 1) * exp(eta_g)"
      )
    )
  ),
  survival_AFT = list(
    outcome_type = "survival",
    stan_func = "lognormal_lpdf",
    extra_params = list(list(name = "sigma", type = "positive")),
    links = list(
      identity = list(
        stan_link = "identity_fn",
        eq = "log(T) | g, X ~ Normal(eta_g, sigma_g)"
      )
    )
  )
)

# ===========================================================================
# Constructor
# ===========================================================================

#' Define a Principal Stratification Model
#'
#' Creates a model specification for principal stratification analysis.
#' No data is required at this stage -- the specification is purely symbolic.
#' Use \code{\link{fit}} to estimate the model with data.
#'
#' @param S.formula formula for the stratum model (e.g., \code{Z + D ~ X1 + X2}).
#' @param Y.formula formula for the outcome model (e.g., \code{Y ~ X1 + X2}).
#' @param Y.family a \code{family} object (e.g., \code{gaussian()}, \code{survival("Cox")}).
#' @param strata strata definition: character vector, list of vectors, or list of lists.
#' @param ER exclusion restriction: character vector of strata names or logical vector.
#' @param prior_intercept,prior_coefficient,prior_sigma,prior_alpha,prior_lambda,prior_theta
#'   prior distributions for model parameters.
#' @param survival.time.points number of time points for survival outcomes.
#'
#' @return An object of class \code{PStrataModel}.
#' @examples
#' # Non-compliance with three strata (never-taker, complier, always-taker)
#' model <- PStrataModel(
#'   S.formula = Z + D ~ 1,
#'   Y.formula = Y ~ 1,
#'   Y.family  = gaussian(),
#'   strata    = c(n = "00", c = "01", a = "11"),
#'   ER        = c("n", "a")
#' )
#' print(model)
#' summary(model)
#'
#' \donttest{
#' # Fit the model (requires rstan and C++ compiler)
#' data(sim_data_normal)
#' ps_fit <- fit(model, data = sim_data_normal, chains = 2, iter = 500)
#' summary(ps_fit)
#' plot(ps_fit)
#'
#' # Extract potential outcomes and contrasts
#' est <- estimate(ps_fit)
#' summary(est)
#' plot(est)
#'
#' ctr <- contrast(ps_fit)
#' summary(ctr)
#' plot(ctr)
#' }
#' @export
PStrataModel <- function(
    S.formula,
    Y.formula,
    Y.family,
    strata,
    ER = NULL,
    prior_intercept = prior_flat(),
    prior_coefficient = prior_normal(),
    prior_sigma = prior_inv_gamma(),
    prior_alpha = prior_inv_gamma(),
    prior_lambda = prior_inv_gamma(),
    prior_theta = prior_normal(),
    survival.time.points = 50
) {
  cl <- match.call()

  # Parse formulas symbolically
  s_parsed <- parse_formula_symbolic(S.formula)
  y_parsed <- parse_formula_symbolic(Y.formula)

  treatment_name <- s_parsed$lhs_names[1]
  postrand_names <- s_parsed$lhs_names[-1]
  outcome_name   <- y_parsed$lhs_names[1]
  event_name     <- if (length(y_parsed$lhs_names) > 1) y_parsed$lhs_names[2] else NULL

  # Parse strata
  strata_parsed <- parse_strata_def(strata)

  # Resolve ER
  er_list <- resolve_ER_model(ER, strata_parsed$strata_names)

  # Cross-validate formula vs strata
  if (length(postrand_names) != strata_parsed$num_postrand) {
    stop(
      "Formula has ", length(postrand_names),
      " post-randomization variable(s) but strata defines ",
      strata_parsed$num_postrand, "."
    )
  }

  # Resolve family
  family_name <- Y.family$family
  family_link <- Y.family$link
  family_info <- resolve_family_info(family_name, family_link)

  is_survival <- family_info$outcome_type == "survival"

  if (is_survival && is.null(event_name)) {
    stop("Survival family requires an event indicator in Y.formula (e.g., Y + delta ~ ...).")
  }

  # Compute outcome groups for display
  outcome_groups <- compute_outcome_groups(
    strata_parsed$strata_names, strata_parsed$strata_def,
    er_list, strata_parsed$num_treatment
  )

  # Collect relevant priors
  priors <- list(intercept = prior_intercept, coefficient = prior_coefficient)
  for (ep in family_info$extra_params) {
    priors[[ep$name]] <- switch(ep$name,
      sigma  = prior_sigma,
      alpha  = prior_alpha,
      lambda = prior_lambda,
      theta  = prior_theta
    )
  }

  res <- list(
    call               = cl,
    stratum_formula     = S.formula,
    outcome_formula     = Y.formula,
    outcome_family      = Y.family,
    treatment_name      = treatment_name,
    postrand_names      = postrand_names,
    outcome_name        = outcome_name,
    event_name          = event_name,
    strata_names        = strata_parsed$strata_names,
    strata_def          = strata_parsed$strata_def,
    exclusion_restriction = er_list,
    num_strata          = length(strata_parsed$strata_names),
    num_treatment       = strata_parsed$num_treatment,
    num_postrand        = length(postrand_names),
    family_info         = family_info,
    priors              = priors,
    outcome_groups      = outcome_groups,
    is_survival         = is_survival,
    survival_time_points = if (is_survival) survival.time.points else NULL,
    s_random_effects    = s_parsed$random_effects,
    y_random_effects    = y_parsed$random_effects,
    s_terms             = s_parsed$rhs_terms,
    y_terms             = y_parsed$rhs_terms
  )

  class(res) <- if (is_survival) {
    c("PStrataModelSurvival", "PStrataModel")
  } else {
    "PStrataModel"
  }

  res
}

# ===========================================================================
# Symbolic formula parsing (no data needed)
# ===========================================================================

#' @noRd
parse_formula_symbolic <- function(formula) {
  lhs <- if (length(formula) == 3) formula[[2]] else NULL
  lhs_names <- extract_lhs_names(lhs)

  bars <- reformulas::findbars(formula)
  random_effects <- NULL

  if (!is.null(bars) && length(bars) > 0) {
    random_effects <- lapply(bars, function(bar) {
      group <- deparse(bar[[3]])
      terms_expr <- bar[[2]]
      terms_str <- deparse(terms_expr)

      tf <- stats::terms(stats::reformulate(terms_str))
      has_intercept <- attr(tf, "intercept") == 1L
      has_slopes <- length(attr(tf, "term.labels")) > 0

      list(
        formula_str  = paste0("(", terms_str, " | ", group, ")"),
        group        = group,
        has_intercept = has_intercept,
        has_slopes   = has_slopes
      )
    })
  }

  # Extract fixed-effect term labels from RHS
  fixed_formula <- reformulas::nobars(formula)
  rhs_terms <- if (length(fixed_formula) == 3) {
    tf <- stats::terms(fixed_formula)
    has_intercept <- attr(tf, "intercept") == 1L
    term_labels <- attr(tf, "term.labels")
    list(has_intercept = has_intercept, covariates = term_labels)
  } else {
    list(has_intercept = FALSE, covariates = character(0))
  }

  list(
    lhs_names      = lhs_names,
    rhs_terms      = rhs_terms,
    random_effects = random_effects
  )
}

#' Recursively extract symbol names from a formula LHS
#' @noRd
extract_lhs_names <- function(node) {
  if (is.null(node)) return(character(0))
  if (is.name(node)) return(as.character(node))
  if (is.call(node)) {
    return(unlist(unique(lapply(node[-1], extract_lhs_names))))
  }
  character(0)
}

# ===========================================================================
# Strata parsing (3 input forms)
# ===========================================================================

#' @noRd
parse_strata_def <- function(strata) {
  if (is.character(strata)) {
    parse_strata_string(strata)
  } else if (is.list(strata)) {
    first <- strata[[1]]
    if (is.list(first)) {
      parse_strata_list_of_lists(strata)
    } else {
      parse_strata_list_of_vectors(strata)
    }
  } else {
    stop("strata must be a character vector or a named list.")
  }
}

#' @noRd
parse_strata_string <- function(strata) {
  strata_clean <- gsub("\\*$", "", strata)
  strata_names <- names(strata)
  if (is.null(strata_names)) strata_names <- strata_clean
  strata_names <- ifelse(strata_names == "", strata_clean, strata_names)

  strata_def <- lapply(strata_clean, function(s) {
    parts <- strsplit(s, "\\|")[[1]]
    lapply(parts, function(p) as.integer(strsplit(p, "")[[1]]))
  })
  names(strata_def) <- strata_names

  list(
    strata_names = strata_names,
    strata_def   = strata_def,
    num_treatment = length(strata_def[[1]][[1]]),
    num_postrand  = length(strata_def[[1]])
  )
}

#' @noRd
parse_strata_list_of_vectors <- function(strata) {
  strata_names <- names(strata)
  if (is.null(strata_names) || any(strata_names == ""))
    stop("All strata must be named.")

  strata_def <- lapply(strata, function(v) list(as.integer(v)))
  names(strata_def) <- strata_names

  list(
    strata_names = strata_names,
    strata_def   = strata_def,
    num_treatment = length(strata[[1]]),
    num_postrand  = 1L
  )
}

#' @noRd
parse_strata_list_of_lists <- function(strata) {
  strata_names <- names(strata)
  if (is.null(strata_names) || any(strata_names == ""))
    stop("All strata must be named.")

  first_names <- names(strata[[1]])
  all_named   <- !is.null(first_names) && all(nzchar(first_names))
  all_unnamed <- is.null(first_names) || all(!nzchar(first_names))

  if (!all_named && !all_unnamed)
    stop("Mixed named and unnamed D variables not supported.")

  if (all_named) {
    for (i in seq_along(strata)) {
      if (!setequal(names(strata[[i]]), first_names))
        stop("All strata must have the same D variable names.")
    }
  }

  strata_def <- lapply(strata, function(s) lapply(s, as.integer))
  names(strata_def) <- strata_names

  list(
    strata_names = strata_names,
    strata_def   = strata_def,
    num_treatment = length(strata[[1]][[1]]),
    num_postrand  = length(strata[[1]])
  )
}

# ===========================================================================
# ER resolution
# ===========================================================================

#' @noRd
resolve_ER_model <- function(ER, strata_names) {
  n <- length(strata_names)
  if (is.null(ER)) return(rep(FALSE, n))
  if (is.logical(ER)) {
    stopifnot(length(ER) == n)
    return(ER)
  }
  if (is.character(ER)) {
    er_list <- rep(FALSE, n)
    for (name in ER) {
      idx <- which(strata_names == name)
      if (length(idx) == 0) warning("ER stratum '", name, "' not found, ignored.")
      else er_list[idx] <- TRUE
    }
    return(er_list)
  }
  stop("ER must be NULL, logical, or character vector.")
}

# ===========================================================================
# Family resolution
# ===========================================================================

#' @noRd
resolve_family_info <- function(family_name, link) {
  reg <- .family_registry[[family_name]]
  if (is.null(reg)) stop("Unsupported family: ", family_name)

  link_info <- reg$links[[link]]
  if (is.null(link_info)) stop("Unsupported link '", link, "' for family '", family_name, "'.")

  list(
    outcome_type = reg$outcome_type,
    stan_func    = reg$stan_func,
    stan_link    = link_info$stan_link,
    extra_params = reg$extra_params,
    equation     = link_info$eq
  )
}

# ===========================================================================
# Outcome group computation
# ===========================================================================

#' @noRd
compute_outcome_groups <- function(strata_names, strata_def, er_list, num_treatment) {
  labels <- character(0)

  for (i in seq_along(strata_names)) {
    if (er_list[i]) {
      labels <- c(labels, strata_names[i])
    } else {
      for (z in seq_len(num_treatment) - 1) {
        labels <- c(labels, paste0(strata_names[i], "(", z, ")"))
      }
    }
  }

  list(
    labels    = labels,
    n_groups  = length(labels),
    er_strata = strata_names[er_list]
  )
}

# ===========================================================================
# print
# ===========================================================================

#' @export
print.PStrataModel <- function(x, ...) {
  cat("PStrataModel\n")
  cat("  Stratum model:", deparse(x$stratum_formula), "\n")
  cat("  Outcome model:", deparse(x$outcome_formula), "\n")
  cat("  Family:", x$outcome_family$family, "(", x$outcome_family$link, ")\n")

  strata_strs <- vapply(seq_along(x$strata_names), function(i) {
    parts <- vapply(x$strata_def[[i]], function(d) paste(d, collapse = ""), character(1))
    paste0(x$strata_names[i], "(", paste(parts, collapse = "|"), ")")
  }, character(1))
  cat("  Strata:", paste(strata_strs, collapse = ", "), "\n")

  er_names <- x$strata_names[x$exclusion_restriction]
  if (length(er_names) > 0)
    cat("  Exclusion restriction:", paste(er_names, collapse = ", "), "\n")

  if (x$is_survival)
    cat("  Time points:", x$survival_time_points, "\n")

  invisible(x)
}

# ===========================================================================
# summary
# ===========================================================================

#' @export
summary.PStrataModel <- function(object, ...) {
  cat("PStrataModel Summary\n")
  cat("====================\n\n")

  # Call
  cat("Call:", deparse(object$call, width.cutoff = 200), "\n\n")

  # Variable roles
  cat("Treatment:", object$treatment_name, "\n")
  cat("Post-randomization:", paste(object$postrand_names, collapse = ", "), "\n")
  if (object$is_survival) {
    cat("Time-to-event outcome:", object$outcome_name, "\n")
    cat("Event indicator:", object$event_name, "(1 = observed, 0 = censored)\n")
  } else {
    cat("Outcome:", object$outcome_name, "\n")
  }

  # Strata table
  cat("\nStrata:\n")
  print_strata_table(object)

  # Stratum model
  cat("\nStratum model (softmax):\n")
  cat("  Stratum:", deparse(object$stratum_formula), "\n")
  if (length(object$s_terms$covariates) > 0)
    cat("  X = [", paste(object$s_terms$covariates, collapse = ", "), "]\n", sep = "")
  re_s <- if (length(object$s_random_effects) > 0) " + r_{s,G}" else ""
  s_eq <- build_linear_predictor("s", object$s_terms, re_s)
  cat("  P(s | X) = softmax(", s_eq, ") for s in {",
      paste(object$strata_names, collapse = ", "), "}\n", sep = "")
  cat("  Priors:\n")
  print_fixed_priors("s", object$s_terms, object$priors)
  print_random_effects(object$s_random_effects, "s")

  # Outcome model
  family_name <- object$outcome_family$family
  family_link <- object$outcome_family$link
  if (object$is_survival) {
    method <- sub("survival_", "", family_name)
    cat("\nOutcome model (survival, ", method, "):\n", sep = "")
  } else {
    cat("\nOutcome model (", family_name, ", ", family_link, " link):\n", sep = "")
  }
  cat("  Outcome:", deparse(object$outcome_formula), "\n")
  if (length(object$y_terms$covariates) > 0)
    cat("  X = [", paste(object$y_terms$covariates, collapse = ", "), "]\n", sep = "")

  re_g <- if (length(object$y_random_effects) > 0) " + r_{g,G}" else ""
  g_eq <- build_linear_predictor("g", object$y_terms, re_g)
  cat("  eta_g = ", g_eq, "\n", sep = "")

  groups_str <- paste(object$outcome_groups$labels, collapse = ", ")
  cat("  ", object$family_info$equation, " for g in {", groups_str, "}\n", sep = "")

  if (length(object$outcome_groups$er_strata) > 0) {
    cat("  Note:", paste(object$outcome_groups$er_strata, collapse = ", "),
        "shared across treatment arms (ER)\n")
  }

  if (object$is_survival)
    cat("  Time points:", object$survival_time_points, "(determined at fit time)\n")

  cat("  Priors:\n")
  print_fixed_priors("g", object$y_terms, object$priors)
  for (ep in object$family_info$extra_params) {
    cat("    ", paste0(ep$name, "_g ~ "), format_prior(object$priors[[ep$name]]), "\n", sep = "")
  }
  print_random_effects(object$y_random_effects, "g")

  # Available estimands
  cat("\nAvailable estimands after fitting:\n")
  if (object$is_survival) {
    cat("  survival probability, RACE (restricted average causal effect)\n")
  } else {
    cat("  mean outcome per stratum and treatment arm\n")
  }

  invisible(object)
}

# ===========================================================================
# Display helpers
# ===========================================================================

#' Build the linear predictor string based on what terms are present
#' @noRd
build_linear_predictor <- function(idx, terms, re_str) {
  parts <- character(0)
  if (terms$has_intercept) parts <- c(parts, paste0("a_", idx))
  if (length(terms$covariates) > 0) parts <- c(parts, paste0("b_", idx, " * X"))
  eq <- paste(parts, collapse = " + ")
  if (nzchar(re_str)) eq <- paste0(eq, re_str)
  eq
}

#' Print prior lines for fixed effects, respecting intercept/covariate presence
#' @noRd
print_fixed_priors <- function(idx, terms, priors) {
  if (terms$has_intercept)
    cat("    a_", idx, ": intercept ~ ", format_prior(priors$intercept), "\n", sep = "")
  if (length(terms$covariates) > 0)
    cat("    b_", idx, ": coefficients ~ ", format_prior(priors$coefficient), "\n", sep = "")
}

#' @noRd
format_prior <- function(prior) {
  if (prior$name == "flat") return("flat")
  params <- paste(unlist(prior$args), collapse = ", ")
  paste0(prior$name, "(", params, ")")
}

#' @noRd
print_strata_table <- function(object) {
  sd   <- object$strata_def
  sn   <- object$strata_names
  pn   <- object$postrand_names
  tn   <- object$treatment_name
  nt   <- object$num_treatment
  er   <- object$exclusion_restriction
  z_levels <- seq_len(nt) - 1

  # Column headers
  col_h <- character(0)
  for (d in pn) {
    for (z in z_levels) col_h <- c(col_h, paste0(d, "(", tn, "=", z, ")"))
  }
  col_h <- c(col_h, "ER")

  # Row data
  rows <- lapply(seq_along(sn), function(i) {
    vals <- unlist(sd[[i]])
    c(as.character(vals), if (er[i]) "*" else "")
  })

  # Column widths
  all_vals <- do.call(rbind, rows)
  cw <- pmax(nchar(col_h), apply(all_vals, 2, function(x) max(nchar(x))))
  nw <- max(nchar(sn))

  # Header
  hdr <- paste0(
    "  ", formatC("", width = nw), "  ",
    paste(mapply(function(h, w) formatC(h, width = w), col_h, cw), collapse = "  ")
  )
  cat(hdr, "\n")

  # Rows
  for (i in seq_along(sn)) {
    r <- paste0(
      "  ", formatC(sn[i], width = nw), "  ",
      paste(mapply(function(v, w) formatC(v, width = w), rows[[i]], cw), collapse = "  ")
    )
    cat(r, "\n")
  }
}

#' @noRd
print_random_effects <- function(re_list, idx_symbol) {
  if (is.null(re_list) || length(re_list) == 0) return(invisible())

  for (re in re_list) {
    cat("  Random effects ", re$formula_str, ":\n", sep = "")
    sym <- paste0("r_{", idx_symbol, ",", re$group, "}")
    if (!re$has_slopes) {
      cat("    ", sym, " ~ normal(0, tau), tau ~ inv_gamma(1, 1)\n", sep = "")
    } else {
      if (re$has_intercept) {
        cat("    ", sym, " has intercept and slope components\n", sep = "")
      } else {
        cat("    ", sym, " has slope components\n", sep = "")
      }
      cat("    each component ~ normal(0, tau), tau ~ inv_gamma(1, 1)\n")
    }
  }
}
