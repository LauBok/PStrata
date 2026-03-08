#' \bold{Stan} Code for \pkg{PStrata} Models
#'
#' Generate the \bold{Stan} code corresponding to the model,
#' which is read by \bold{Stan} to do sampling.
#'
#' @param PSobject an object of class \code{PSobject}
#' @param filename (optional) string. If not \code{NULL}, the stan file will be saved via
#' \code{\link{cat}} in a text file named after the string supplied.
#' @param debug only for testing in development mode. Will be removed in future release.
#'
#' @return A string, which can be printed on screen using \code{\link{cat}}.
#'
#' @export
make_stancode <- function(PSobject, filename = NULL, debug = FALSE) {
  ctx <- build_stan_context(PSobject, debug)

  blocks <- c(
    stan_comment_lines(),
    " ",
    stan_functions_block(ctx),
    " ",
    stan_data_block(ctx),
    " ",
    stan_transformed_data_block(ctx),
    " ",
    stan_parameters_block(ctx),
    " ",
    stan_transformed_parameters_block(ctx),
    " ",
    stan_model_block(ctx),
    " ",
    stan_generated_quantities_block(ctx)
  )

  full_text <- paste0(blocks, collapse = "\n")

  if (!is.null(filename)) {
    if (!stringr::str_ends(filename, ".stan"))
      filename <- paste0(filename, ".stan")
    cat(full_text, file = filename)
  }

  full_text
}

# =============================================================================
# Context builder
# =============================================================================

#' Read an inst/ config file, handling debug vs installed paths
#' @noRd
read_inst_file <- function(filename, debug) {
  if (debug) {
    return(readLines(file.path("inst", filename)))
  }
  pkg_path <- path.package("PStrata")
  # Installed packages have files directly under pkg_path;

  # devtools::load_all() keeps them under inst/
  path <- file.path(pkg_path, filename)
  if (!file.exists(path)) {
    path <- file.path(pkg_path, "inst", filename)
  }
  readLines(path)
}

#' Parse family_info.txt to extract outcome type, Stan function, and extra params
#' @noRd
parse_family_info <- function(lines, Y_family) {
  Y_type <- NULL
  func_family <- NULL
  parameter_list <- list()

  for (line in lines) {
    parts <- stringr::str_split(line, " +")[[1]]
    if (parts[1] != Y_family) next

    Y_type <- parts[2]
    func_family <- parts[3]
    n_params <- as.integer(parts[4])
    if (n_params > 0) {
      for (i in seq_len(n_params)) {
        parameter_list[[i]] <- list(
          name = parts[2 * i + 3],
          type = parts[2 * i + 4]
        )
      }
    }
    break
  }

  list(Y_type = Y_type, func_family = func_family, parameter_list = parameter_list)
}

#' Parse link_info.txt to find the Stan inverse-link function name
#' @noRd
parse_link_info <- function(lines, Y_family, Y_link) {
  for (line in lines) {
    parts <- stringr::str_split(line, " +")[[1]]
    if (parts[1] == Y_family && parts[2] == Y_link) return(parts[3])
  }
  NULL
}

#' Build the prior list from PSobject for Stan code generation
#' @noRd
build_prior_list <- function(PSobject) {
  prior_names <- c("intercept", "coefficient", "sigma", "alpha", "lambda", "theta")
  lapply(prior_names, function(name) {
    prior_curr <- PSobject[[paste0("prior_", name)]]
    list(
      type  = name,
      func  = prior_curr$name,
      param = unlist(prior_curr$args)
    )
  })
}

#' Assemble all derived quantities needed for Stan code generation
#' @noRd
build_stan_context <- function(PSobject, debug) {
  SZDG_max <- apply(PSobject$SZDG_table, 2, max)
  Y_family <- PSobject$Y.family$family
  Y_link   <- PSobject$Y.family$link

  family_info <- parse_family_info(
    read_inst_file("family_info.txt", debug), Y_family
  )
  func_link <- parse_link_info(
    read_inst_file("link_info.txt", debug), Y_family, Y_link
  )

  list(
    PSobject       = PSobject,
    SZDG_max       = SZDG_max,
    S_re           = length(PSobject$S.formula$random_eff_list),
    Y_re           = length(PSobject$Y.formula$random_eff_list),
    Y_family       = Y_family,
    Y_link         = Y_link,
    Y_type         = family_info$Y_type,
    func_family    = family_info$func_family,
    func_link      = func_link,
    parameter_list = family_info$parameter_list,
    prior_list     = build_prior_list(PSobject),
    func_imp_lines = read_inst_file("function_implement.txt", debug)
  )
}

# =============================================================================
# Random effect line generators (shared across multiple blocks)
# =============================================================================

#' Generate Stan data declarations for random effects
#' @param prefix "S" for stratum model, "G" for outcome model
#' @param re_count Number of random effect groups
#' @noRd
stan_re_data_lines <- function(prefix, re_count) {
  if (re_count == 0) return(character(0))
  model_desc <- if (prefix == "S") "principal stratum" else "outcome"
  lines <- paste0("    // random effect for ", model_desc, " model")
  for (i in seq_len(re_count)) {
    p <- paste0("P", prefix, "_RE_", i)
    n <- paste0("N", prefix, "_RE_", i)
    x <- paste0("X", prefix, "_RE_", i)
    lines <- c(lines,
      paste0("    int<lower=0> ", p, "; // number of random effect terms"),
      paste0("    int<lower=0> ", n, "; // number of levels"),
      paste0("    matrix[N, ", p, "*", n, "] ", x, "; // model matrix for random effect")
    )
  }
  lines
}

#' Generate Stan parameter declarations for random effects
#' @noRd
stan_re_param_lines <- function(prefix, max_dim, re_count) {
  if (re_count == 0) return(character(0))
  model_desc <- if (prefix == "S") "principal stratum" else "outcome"
  lines <- paste0("    // random effect for ", model_desc, " model")
  for (i in seq_len(re_count)) {
    lines <- c(lines,
      paste0("    matrix[", max_dim, ", P", prefix, "_RE_", i,
             "*N", prefix, "_RE_", i, "] beta_", prefix, "_RE_", i, ";"),
      paste0("    real<lower=0> tau_", prefix, "_RE_", i,
             "[", max_dim, ", P", prefix, "_RE_", i, "];")
    )
  }
  lines
}

#' Generate transformed parameter declarations for random effects
#' @noRd
stan_re_transformed_decl <- function(prefix, max_dim, re_count) {
  if (re_count == 0) return(character(0))
  lines <- character(0)
  for (i in seq_len(re_count)) {
    lines <- c(lines,
      paste0("    matrix[N", prefix, "_RE_", i, ", P", prefix, "_RE_", i,
             "] M_beta_", prefix, "_RE_", i, "[", max_dim, "];")
    )
  }
  lines
}

#' Generate transformed parameter body (reshape) for random effects
#' @param prefix "S" or "G"
#' @param re_count Number of RE groups
#' @param loop_bound The Stan for-loop upper bound
#' @param body_count Number of assignment lines to generate inside the loop
#' @noRd
stan_re_transformed_body <- function(prefix, re_count, loop_bound, body_count) {
  if (re_count == 0) return(character(0))
  lines <- character(0)
  for (i in seq_len(re_count)) {
    lines <- c(lines, paste0("    for (i in 1:", loop_bound, "){"))
    for (j in seq_len(body_count)) {
      lines <- c(lines,
        paste0("        M_beta_", prefix, "_RE_", i, "[", j,
               "] = to_matrix(beta_", prefix, "_RE_", i, "[", j,
               "], N", prefix, "_RE_", i, ", P", prefix, "_RE_", i, ", 0);")
      )
    }
    lines <- c(lines, "    }")
  }
  lines
}

#' Generate model-block priors for random effects
#' @noRd
stan_re_prior_lines <- function(prefix, max_dim, re_count) {
  if (re_count == 0) return(character(0))
  lines <- character(0)
  for (i in seq_len(re_count)) {
    lines <- c(lines,
      paste0("    for (i in 1:", max_dim, ") {"),
      paste0("        tau_", prefix, "_RE_", i, "[i] ~ inv_gamma(1, 1);"),
      paste0("        for (j in 1:P", prefix, "_RE_", i, ") {"),
      paste0("            M_beta_", prefix, "_RE_", i,
             "[i][:, j] ~ normal(0, tau_", prefix, "_RE_", i, "[i, j]);"),
      "        }",
      "    }"
    )
  }
  lines
}

#' Generate the RE linear predictor addition terms (e.g. " + XS_RE_1[n] * beta_S_RE_1[s-1]'")
#' @noRd
stan_re_linear_terms <- function(prefix, re_count, row_idx, coef_idx) {
  if (re_count == 0) return("")
  paste(vapply(seq_len(re_count), function(i) {
    paste0(" + X", prefix, "_RE_", i, "[", row_idx, "] * beta_", prefix, "_RE_", i, "[", coef_idx, "]'")
  }, character(1)), collapse = "")
}

#' Generate RE matrix linear terms (no row index, for generated quantities)
#' @noRd
stan_re_matrix_terms <- function(prefix, re_count) {
  if (re_count == 0) return("")
  paste(vapply(seq_len(re_count), function(i) {
    paste0(" + X", prefix, "_RE_", i, " * beta_", prefix, "_RE_", i, "'")
  }, character(1)), collapse = "")
}

# =============================================================================
# Stan block generators
# =============================================================================

#' @noRd
stan_comment_lines <- function() {
  c(
    "// Automatically generated by PStrata 0.0.4",
    paste0("// at ", Sys.time(), " ", Sys.timezone())
  )
}

#' @noRd
stan_functions_block <- function(ctx) {
  lines <- character(0)
  aim <- FALSE
  for (line in ctx$func_imp_lines) {
    if (line == "<<<>>>") aim <- FALSE
    if (!aim) {
      parts <- stringr::str_split(line, " ")[[1]]
      if (parts[1] == "<<<") {
        if (parts[2] %in% c(ctx$func_family, ctx$func_link))
          aim <- TRUE
      }
      next
    }
    lines <- c(lines, paste0("    ", line))
  }
  if (length(lines) > 0) c("functions {", lines, "}") else character(0)
}

#' @noRd
stan_data_block <- function(ctx) {
  mx <- ctx$SZDG_max
  lines <- c(
    "data {",
    "    int<lower=0> N; // number of observations",
    "    int<lower=0> PS; // number of predictors for principal stratum model",
    "    int<lower=0> PG; // number of predictors for outcome model",
    paste0("    int<lower=0, upper=", mx[2], "> Z[N]; // treatment arm"),
    paste0("    int<lower=0, upper=", mx[3], "> D[N]; // post randomization confounding variable")
  )

  Y_lines <- switch(ctx$Y_type,
    "real"     = "    real Y[N]; // real outcome",
    "positive" = "    real<lower=0> Y[N]; // positive outcome",
    "binary"   = "    int<lower=0, upper=1> Y[N]; // binary outcome",
    "count"    = "    int<lower=0> Y[N]; // count outcome",
    "survival" = c(
      "    real<lower=0> Y[N]; // survival outcome",
      "    int<lower=0, upper=1> delta[N]; // event indicator",
      "    int<lower=0> T; // number of time points for evaluation",
      "    real<lower=0> time[T]; // time points"
    )
  )

  c(lines, Y_lines,
    "    matrix[N, PS] XS; // model matrix for principal stratum model",
    "    matrix[N, PG] XG; // model matrix for outcome model",
    stan_re_data_lines("S", ctx$S_re),
    stan_re_data_lines("G", ctx$Y_re),
    "}"
  )
}

#' @noRd
stan_transformed_data_block <- function(ctx) {
  mx <- ctx$SZDG_max
  SG_table <- unique(ctx$PSobject$SZDG_table[, c("S", "G")])

  lines <- c(
    "transformed data {",
    paste0("    int S[", mx[4], "];")
  )
  for (i in seq_len(nrow(SG_table))) {
    lines <- c(lines,
      paste0("   S[", SG_table[i, "G"], "] = ", SG_table[i, "S"] + 1, ";")
    )
  }
  c(lines, "}")
}

#' @noRd
stan_parameters_block <- function(ctx) {
  mx <- ctx$SZDG_max
  lines <- c(
    "parameters {",
    paste0("    matrix[", mx['S'], ", PS] beta_S; // coefficients for principal stratum model"),
    paste0("    matrix[", mx['G'], ", PG] beta_G; // coefficients for outcome model"),
    stan_re_param_lines("S", mx['S'], ctx$S_re),
    stan_re_param_lines("G", mx['G'], ctx$Y_re)
  )

  for (param in ctx$parameter_list) {
    constraint <- if (param$type == "positive") "<lower=0> " else " "
    lines <- c(lines,
      paste0("    real", constraint, param$name, "[", mx['G'], "];")
    )
  }

  c(lines, "}")
}

#' @noRd
stan_transformed_parameters_block <- function(ctx) {
  mx <- ctx$SZDG_max
  c(
    "transformed parameters {",
    stan_re_transformed_decl("S", mx['S'], ctx$S_re),
    stan_re_transformed_decl("G", mx['G'], ctx$Y_re),
    stan_re_transformed_body("S", ctx$S_re, mx['S'], mx['S']),
    stan_re_transformed_body("G", ctx$Y_re, mx['S'], mx['G']),
    "}"
  )
}

#' @noRd
stan_model_block <- function(ctx) {
  mx <- ctx$SZDG_max
  PSobject <- ctx$PSobject

  lines <- c(
    "model {",
    "    // random effect",
    stan_re_prior_lines("S", mx['S'], ctx$S_re),
    stan_re_prior_lines("G", mx['G'], ctx$Y_re),
    "    // prior",
    stan_prior_lines(ctx),
    "    // model",
    stan_likelihood_lines(ctx)
  )

  c(lines, "}")
}

#' Generate prior statements for the model block
#' @noRd
stan_prior_lines <- function(ctx) {
  lines <- character(0)
  for (prior in ctx$prior_list) {
    if (prior$func == "flat") next
    params_str <- paste(prior$param, collapse = ", ")

    if (prior$type == "intercept") {
      lines <- c(lines,
        paste0("    beta_S[:, 1] ~ ", prior$func, "(", params_str, ");"),
        paste0("    beta_G[:, 1] ~ ", prior$func, "(", params_str, ");")
      )
    } else if (prior$type == "coefficient") {
      lines <- c(lines,
        "    if (PS >= 2)",
        paste0("        to_vector(beta_S[:, 2:PS]) ~ ", prior$func, "(", params_str, ");"),
        "    if (PG >= 2)",
        paste0("        to_vector(beta_G[:, 2:PG]) ~ ", prior$func, "(", params_str, ");")
      )
    } else {
      is_model_param <- any(vapply(ctx$parameter_list,
        function(p) p$name == prior$type, logical(1)))
      if (is_model_param) {
        lines <- c(lines,
          paste0("    ", prior$type, " ~ ", prior$func, "(", params_str, ");")
        )
      }
    }
  }
  lines
}

#' Generate the main likelihood loop for the model block
#' @noRd
stan_likelihood_lines <- function(ctx) {
  mx <- ctx$SZDG_max
  PSobject <- ctx$PSobject
  S_re_term <- stan_re_linear_terms("S", ctx$S_re, "n", "s-1")

  lines <- c(
    "    for (n in 1:N) {",
    "        int length;",
    paste0("        real log_prob[", mx['S'] + 1, "];"),
    "        log_prob[1] = 0;",
    paste0("        for (s in 2:", mx['S'] + 1, ") {"),
    paste0("            log_prob[s] = XS[n] * beta_S[s-1]'", S_re_term, ";"),
    "        }"
  )

  ZD_comb <- unique(PSobject$SZDG_table[, c("Z", "D")])

  # length assignment block
  b_else <- FALSE
  for (r in seq_len(nrow(ZD_comb))) {
    z <- ZD_comb[r, "Z"]; d <- ZD_comb[r, "D"]
    n_groups <- sum(PSobject$SZDG_table[, "Z"] == z & PSobject$SZDG_table[, "D"] == d)
    lines <- c(lines,
      paste0("        ", if (b_else) "else ", "if (Z[n] == ", z, " && D[n] == ", d, ")"),
      paste0("            length = ", n_groups, ";")
    )
    b_else <- TRUE
  }

  lines <- c(lines,
    "        {",
    "            real log_l[length];"
  )

  # likelihood computation block
  b_else <- FALSE
  for (r in seq_len(nrow(ZD_comb))) {
    z <- ZD_comb[r, "Z"]; d <- ZD_comb[r, "D"]
    SG_table <- PSobject$SZDG_table[
      PSobject$SZDG_table[, "Z"] == z & PSobject$SZDG_table[, "D"] == d, , drop = FALSE
    ]
    lines <- c(lines,
      paste0("            ", if (b_else) "else ", "if (Z[n] == ", z, " && D[n] == ", d, ") {"),
      paste0("                // Z:", z, " D:", d, " S:", paste(SG_table[, "S"], collapse = "/"))
    )

    for (sg in seq_len(nrow(SG_table))) {
      s <- SG_table[sg, "S"]
      g <- SG_table[sg, "G"]
      Y_re_term <- stan_re_linear_terms("G", ctx$Y_re, "n", g)

      str_param <- ""
      if (length(ctx$parameter_list) > 0) {
        param_refs <- paste(vapply(ctx$parameter_list,
          function(p) paste0(p$name, "[", g, "]"), character(1)), collapse = ", ")
        str_param <- paste0(", ", param_refs,
          if (ctx$Y_type == "survival") ", delta[n]" else "")
      }

      lines <- c(lines,
        paste0("                log_l[", sg, "] = log_prob[", s + 1, "] + ",
               ctx$func_family, "(Y[n] | ",
               ctx$func_link, "(XG[n] * beta_G[", g, "]'", Y_re_term, ")",
               str_param, ");")
      )
    }
    lines <- c(lines, "            }")
    b_else <- TRUE
  }

  c(lines,
    "            target += log_sum_exp(log_l) - log_sum_exp(log_prob);",
    "        }",
    "    }"
  )
}

# =============================================================================
# Generated quantities block
# =============================================================================

#' @noRd
stan_generated_quantities_block <- function(ctx) {
  mx <- ctx$SZDG_max

  if (ctx$Y_type == "survival") {
    stan_gq_survival(ctx, mx)
  } else {
    stan_gq_ordinary(ctx, mx)
  }
}

#' Generated quantities for non-survival outcomes
#' @noRd
stan_gq_ordinary <- function(ctx, mx) {
  S_re_mat <- stan_re_matrix_terms("S", ctx$S_re)
  Y_re_term_fn <- function(idx_var) stan_re_linear_terms("G", ctx$Y_re, "i", idx_var)

  lines <- c(
    "generated quantities {",
    paste0("    vector[", mx['S'] + 1, "] strata_prob; // the probability of being in each stratum"),
    paste0("    vector[", mx['G'], "] mean_effect; // mean response"),
    "    {",
    paste0("        matrix[N, ", mx['S'] + 1, "] log_prob;"),
    paste0("        vector[", mx['G'], "] numer;"),
    paste0("        matrix[N, ", mx['G'], "] expected_mean;"),
    "        for (i in 1:N)",
    paste0("            for (j in 1:", mx['G'], ")"),
    paste0("                expected_mean[i, j] = ", ctx$func_link,
           "(XG[i] * beta_G[j]'", Y_re_term_fn("j"), ");"),
    "        log_prob[:, 1] = rep_vector(0, N);",
    paste0("        log_prob[:, 2:", mx['S'] + 1, "] = XS * beta_S'", S_re_mat, ";"),
    "        for (n in 1:N) {",
    "            log_prob[n] -= log_sum_exp(log_prob[n]);",
    "        }",
    paste0("        for (s in 1:", mx['S'] + 1, ") strata_prob[s] = mean(exp(log_prob[:, s]));"),
    paste0("        for (g in 1:", mx['G'], ") {"),
    "            numer[g] = mean(expected_mean[:, g] .* exp(log_prob[:, S[g]]));",
    "            mean_effect[g] = numer[g] / strata_prob[S[g]];",
    "        }",
    "    }",
    "}"
  )
  lines
}

#' Generated quantities for survival outcomes
#' @noRd
stan_gq_survival <- function(ctx, mx) {
  S_re_mat <- stan_re_matrix_terms("S", ctx$S_re)
  Y_re_term_fn <- function(idx_var) stan_re_linear_terms("G", ctx$Y_re, "i", idx_var)

  lines <- c(
    "generated quantities {",
    paste0("    vector[", mx['S'] + 1, "] strata_prob; // the probability of being in each stratum"),
    paste0("    matrix[", mx['G'], ", T] mean_surv_prob; // mean survival probability"),
    paste0("    matrix[", mx['G'], ", T] mean_RACE; // mean restricted average causal effect"),
    "    {",
    paste0("        matrix[N, ", mx['S'] + 1, "] log_prob;"),
    paste0("        matrix[", mx['G'], ", T] numer_surv_prob;"),
    paste0("        matrix[", mx['G'], ", T] numer_RACE;"),
    paste0("        matrix[N, T] expected_surv_prob[", mx['G'], "];"),
    paste0("        matrix[N, T] expected_RACE[", mx['G'], "];"),
    "        for (i in 1:N)",
    paste0("            for (j in 1:", mx['G'], ")"),
    "                for (t in 1:T) {",
    paste0("                    real mu = ", ctx$func_link,
           "(XG[i] * beta_G[j]'", Y_re_term_fn("j"), ");")
  )

  if (ctx$Y_family == "survival_Cox") {
    lines <- c(lines,
      "                    expected_surv_prob[j][i,t] = exp(-exp(mu) * pow(time[t], exp(theta[j])));",
      "                    expected_RACE[j][i,t] = exp(-theta[j] - mu * exp(-theta[j])) * tgamma(exp(-theta[j])) * gamma_p(exp(-theta[j]), exp(mu) * pow(time[t], exp(theta[j])));"
    )
  } else if (ctx$Y_family == "survival_AFT") {
    lines <- c(lines,
      "                    expected_surv_prob[j][i,t] = 1 - normal_cdf(log(time[t]), mu, sigma[j]);",
      "                    expected_RACE[j][i,t] = time[t] - normal_cdf(log(time[t]), mu, sigma[j]) + exp(mu + sigma[j] * sigma[j] / 2) * normal_cdf(log(time[t]), mu + sigma[j] * sigma[j], sigma[j]);"
    )
  }

  c(lines,
    "                }",
    "        log_prob[:, 1] = rep_vector(0, N);",
    paste0("        log_prob[:, 2:", mx['S'] + 1, "] = XS * beta_S'", S_re_mat, ";"),
    "        for (n in 1:N) {",
    "            log_prob[n] -= log_sum_exp(log_prob[n]);",
    "        }",
    paste0("        for (s in 1:", mx['S'] + 1, ") strata_prob[s] = mean(exp(log_prob[:, s]));"),
    paste0("        for (g in 1:", mx['G'], ") {"),
    "            for (t in 1:T) {",
    "                numer_surv_prob[g, t] = mean(expected_surv_prob[g][:, t] .* exp(log_prob[:, S[g]]));",
    "                numer_RACE[g, t] = mean(expected_RACE[g][:, t] .* exp(log_prob[:, S[g]]));",
    "                mean_surv_prob[g, t] = numer_surv_prob[g, t] / strata_prob[S[g]];",
    "                mean_RACE[g, t] = numer_RACE[g, t] / strata_prob[S[g]];",
    "            }",
    "        }",
    "    }",
    "}"
  )
}
