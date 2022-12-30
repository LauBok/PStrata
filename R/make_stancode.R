#' \bold{Stan} Code for \pkg{PStrata} Models
#' 
#' Generate the \bold{Stan} code corresponding to the model,
#' which is read by \bold{Stan} to do sampling.
#' 
#' @param PSobject an object of class \code{PSobject}
#' @param filename (optional) string. If not \code{NULL}, the stan file will be saved via
#' \code{\link{cat}} in a text file named after the string supplied.
#'
#' @return A string, which can be printed on screen using \code{\link{cat}}.
#'
#' @export
make_stancode <- function(PSobject, filename = NULL) {
  ## initialization ----------
  SZDG_max <- apply(PSobject$SZDG_table, 2, max)
  S_re <- length(PSobject$S.formula$random_eff_list)
  Y_re <- length(PSobject$Y.formula$random_eff_list)
  Y_family <- PSobject$Y.family$family
  Y_link <- PSobject$Y.family$link
  Y_type <- NULL
  func_family <- NULL
  func_link <- NULL
  parameter_list <- list()
  prior_list <- list()
  
  ### prior ----------
  prior_names <- c("intercept", "coefficient", "sigma", "alpha", "lambda", "theta")
  for(name in prior_names){
    prior_curr <- eval(parse(text = paste0("PSobject$prior_", name)))
    prior_args <- prior_curr$args
    prior_list <- append(prior_list, list(list(
      type = name,
      func = prior_curr$name,
      param = unlist(prior_args)
    )))
  }
  
  ### family info ----------
  family_info_lines <- readLines("inst/family_info.txt")
  splitted_family_info_lines <- stringr::str_split(family_info_lines, " +")
  for (splitted_line in splitted_family_info_lines) {
    if (Y_family == splitted_line[1]) {
      if (as.integer(splitted_line[4]) > 0){
        for (i in 1:as.integer(splitted_line[4])){
          parameter_list <- append(parameter_list, list(list(
            name = splitted_line[2 * i + 3],
            type = splitted_line[2 * i + 4]
          )))
        }
      }
      func_family <- splitted_line[3]
      Y_type <- splitted_line[2]
    }
  }
  
  ### link info ---------
  link_info_lines <- readLines("inst/link_info.txt")
  splitted_link_info_lines <- stringr::str_split(link_info_lines, " +")
  for (splitted_line in splitted_link_info_lines) {
    if (Y_family == splitted_line[1] && Y_link == splitted_line[2]) {
      func_link <- splitted_line[3]
    }
  }
  
  ## functions ----------
  func_imp_lines <- readLines("inst/function_implement.txt")
  func_output_lines <- c()
  aim <- F
  for (line in func_imp_lines) {
    if (line == "<<<>>>") aim <- F
    if (!aim) {
      splitted <- stringr::str_split(line, " ")[[1]]
      if (splitted[1] == "<<<") {
        func_name <- splitted[2]
        if (func_name == func_family || func_name == func_link)
          aim <- T
      }
      next
    }
    else {
      func_output_lines <- c(func_output_lines, paste0("    ", line)) 
    }
  }
  if (length(func_output_lines) > 0)
    func_output_lines <- c("functions {", func_output_lines, "}")
  
  ## stan_data -----------
  stan_data_output_lines <- c(
    "data {",
    "    int<lower=0> N; // number of observations",
    "    int<lower=0> PS; // number of predictors for principal stratum model",
    "    int<lower=0> PG; // number of predictors for outcome model",
    paste0("    int<lower=0, upper=", SZDG_max[2], "> Z[N]; // treatment arm"),
    paste0("    int<lower=0, upper=", SZDG_max[3], "> D[N]; // post randomization confounding variable")
  )
  if (Y_type == "real") {
    Y_line <- "    real Y[N]; // real outcome"
  } else if (Y_type == "positive") {
    Y_line <- "    real<lower=0> Y[N]; // positive outcome"
  } else if (Y_type == "binary") {
    Y_line <- "    int<lower=0, upper=1> Y[N]; // binary outcome"
  } else if (Y_type == "count") {
    Y_line <- "    int<lower=0> Y[N]; // count outcome"
  } else if (Y_type == "survival") {
    Y_line <- c(
      "    real<lower=0> Y[N]; // survival outcome",
      "    int<lower=0, upper=1> delta[N]; // event indicator",
      "    int<lower=0> T; // number of time points for evaluation",
      "    real<lower=0> time[T]; // time points"
    )
  }
  stan_data_output_lines <- c(stan_data_output_lines, Y_line)
  stan_data_output_lines <- c(
    stan_data_output_lines,
    "    matrix[N, PS] XS; // model matrix for principal stratum model",
    "    matrix[N, PG] XG; // model matrix for outcome model"
  )
  if (S_re > 0) {
    stan_data_output_lines <- c(
      stan_data_output_lines,
      "    // random effect for principal stratum model"
    )
    for (i in 1:S_re) {
      tmp_P <- paste0("PS_RE_", i)
      tmp_N <- paste0("NS_RE_", i)
      tmp_X <- paste0("XS_RE_", i)
      stan_data_output_lines <- c(
        stan_data_output_lines, 
        paste0("    int<lower=0> ", tmp_P, "; // number of random effect terms"),
        paste0("    int<lower=0> ", tmp_N, "; // number of levels"),
        paste0("    matrix[N, ", tmp_P, "*", tmp_N, "] " + tmp_X + 
                 "; // model matrix for random effect")
      )
    }
  }
  if (Y_re > 0) {
    stan_data_output_lines <- c(
      stan_data_output_lines,
      "    // random effect for outcome model"
    )
    for (i in 1:Y_re) {
      tmp_P <- paste0("PG_RE_", i)
      tmp_N <- paste0("NG_RE_", i)
      tmp_X <- paste0("XG_RE_", i)
      stan_data_output_lines <- c(
        stan_data_output_lines, 
        paste0("    int<lower=0> ", tmp_P, "; // number of random effect terms"),
        paste0("    int<lower=0> ", tmp_N, "; // number of levels"),
        paste0("    matrix[N, ", tmp_P, "*", tmp_N, "] " + tmp_X + 
                 "; // model matrix for random effect")
      )
    }
  }
  stan_data_output_lines <- c(
    stan_data_output_lines,
    "}"
  )
  
  ## transformed data -------------
  transformed_data_output_lines <- c(
    "transformed data {",
    paste0("    int S[", SZDG_max[4], "];")
  )
  SG_table <- unique(PSobject$SZDG_table[, c("S", "G")])
  for (i in 1:nrow(SG_table)) {
    transformed_data_output_lines <- c(
      transformed_data_output_lines,
      paste0("   S[", SG_table[i, "G"], "] = ", SG_table[i, "S"] + 1, ";")
    )
  }
  transformed_data_output_lines <- c(
    transformed_data_output_lines,
    "}"
  )
  
  ## parameters ---------------
  parameters_output_lines <- c(
    "parameters {",
    paste0("    matrix[", SZDG_max['S'], ", PS] beta_S; // coefficients for principal stratum model"),
    paste0("    matrix[", SZDG_max['G'], ", PG] beta_G; // coefficients for outcome model")
  )
  if (S_re > 0) {
    parameters_output_lines <- c(
      parameters_output_lines,
      "    // random effect for principal stratum model"
    )
    for (i in 1:S_re) {
      parameters_output_lines <- c(
        parameters_output_lines,
        paste0("    matrix[", SZDG_max['S'], ", PS_RE_", i, "*NS_RE_", i, "] beta_S_RE_", i, ";"),
        paste0("    real<lower=0> tau_S_RE_", i, "[", SZDG_max['S'], ", PS_RE_", i, "];")
      )
    }
  }
  if (Y_re > 0) {
    parameters_output_lines <- c(
      parameters_output_lines,
      "    // random effect for outcome model"
    )
    for (i in 1:Y_re) {
      parameters_output_lines <- c(
        parameters_output_lines,
        paste0("    matrix[", SZDG_max['G'], ", PG_RE_", i, "*NG_RE_", i, "] beta_G_RE_", i, ";"),
        paste0("    real<lower=0> tau_G_RE_", i, "[", SZDG_max['G'], ", PG_RE_", i, "];")
      )
    }
  }
  for (parameter in parameter_list) {
    parameters_output_lines <- c(
      parameters_output_lines,
      paste0("    real", ifelse(parameter$type == "positive", "<lower=0> ", " "), 
             parameter$name, "[", SZDG_max['G'], "];")
    )
  }
  parameters_output_lines <- c(
    parameters_output_lines,
    "}"
  )
  
  ## transformed parameters ------------
  transformed_parameters_output_lines <- c(
    "transformed parameters {"
  )
  if (S_re > 0){
    for (i in 1:S_re) {
      transformed_parameters_output_lines <- c(
        transformed_parameters_output_lines,
        paste0("    matrix[NS_RE_", i, ", PS_RE_", i, "] M_beta_S_RE_", i, "[", SZDG_max['S'], "];")
      )
    }
  }
  if (Y_re > 0){
    for (i in 1:Y_re) {
      transformed_parameters_output_lines <- c(
        transformed_parameters_output_lines,
        paste0("    matrix[NG_RE_", i, ", PG_RE_", i, "] M_beta_G_RE_", i, "[", SZDG_max['G'], "];")
      )
    }
  }
  if (S_re > 0){
    for (i in 1:S_re) {
      transformed_parameters_output_lines <- c(
        transformed_parameters_output_lines,
        paste0("    for (i in 1:", SZDG_max['S'], "){")
      )
      for (j in 1:SZDG_max['S']){
        transformed_parameters_output_lines <- c(
          transformed_parameters_output_lines,
          paste0("        M_beta_S_RE_", i, "[", j, "] = to_matrix(beta_S_RE_", i, "[", j, "], NS_RE_", i, 
                 ", PS_RE_", i, ", 0);")
        )
      }
      transformed_parameters_output_lines <- c(
        transformed_parameters_output_lines,
        "    }"
      )
    }
  }
  if (Y_re > 0){
    for (i in 1:Y_re) {
      transformed_parameters_output_lines <- c(
        transformed_parameters_output_lines,
        paste0("    for (i in 1:", SZDG_max['S'], "){")
      )
      for (j in 1:SZDG_max['G']){
        transformed_parameters_output_lines <- c(
          transformed_parameters_output_lines,
          paste0("        M_beta_G_RE_", i, "[", j, "] = to_matrix(beta_G_RE_", i, "[", j, "], NG_RE_", i, 
                 ", PG_RE_", i, ", 0);")
        )
      }
      transformed_parameters_output_lines <- c(
        transformed_parameters_output_lines,
        "    }"
      )
    }
  }
  transformed_parameters_output_lines <- c(
    transformed_parameters_output_lines,
    "}"
  )
  
  ## model --------------------------
  model_output_lines <- c(
    "model {",
    "    // random effect"
  )
  if (S_re > 0) {
    for (i in 1:S_re) {
      model_output_lines <- c(
        model_output_lines,
        paste0("    for (i in 1:", SZDG_max['S'], ") {"),
        paste0("        tau_S_RE_", i, "[i] ~ inv_gamma(1, 1);"),
        paste0("        for (j in 1:PS_RE_", i, ") {"),
        paste0("            M_beta_S_RE_", i, "[i][:, j] ~ normal(0, tau_S_RE_", i, "[i, j]);"),
        "        }",
        "    }"
      )
    }
  }
  if (Y_re > 0) {
    for (i in 1:Y_re) {
      model_output_lines <- c(
        model_output_lines,
        paste0("    for (i in 1:", SZDG_max['G'], ") {"),
        paste0("        tau_G_RE_", i, "[i] ~ inv_gamma(1, 1);"),
        paste0("        for (j in 1:PG_RE_", i, ") {"),
        paste0("            M_beta_G_RE_", i, "[i][:, j] ~ normal(0, tau_G_RE_", i, "[i, j]);"),
        "        }",
        "    }"
      )
    }
  }
  model_output_lines <- c(
    model_output_lines,
    "    // prior"
  )
  for (prior in prior_list) {
    type <- prior$type
    func <- prior$func
    params <- prior$param
    params_str <- paste(params, collapse = ", ")
    if (func == "flat") next
    if (type == "intercept") {
      model_output_lines <- c(
        model_output_lines,
        paste0("    beta_S[:, 1] ~ ", func, "(", params_str, ");"),
        paste0("    beta_G[:, 1] ~ ", func, "(", params_str, ");")
      )
    } else if (type == "coefficient") {
      model_output_lines <- c(
        model_output_lines,
        paste0("    to_vector(beta_S[:, 2:PS]) ~ ", func, "(", params_str, ");"),
        paste0("    to_vector(beta_G[:, 2:PS]) ~ ", func, "(", params_str, ");")
      )
    } else {
      found <- F
      for (parameter in parameter_list)
        if (type == parameter$name) {
          found <- T
          break
        }
      if (found) {
        model_output_lines <- c(
          model_output_lines,
          paste0("    ", type, " ~ ", func, "(", params_str, ");")
        )
      }
    }
  }
  model_output_lines <- c(
    model_output_lines,
    "    // model",
    "    for (n in 1:N) {",
    "        int length;",
    paste0("        real log_prob[", SZDG_max['S'] + 1, "];"),
    "        log_prob[1] = 0;",
    paste0("        for (s in 2:", SZDG_max['S'] + 1, ") {"),
    paste0(
      "            log_prob[s] = XS[n] * beta_S[s-1]'",
      ifelse(S_re == 0, "", paste(
        sapply(1:S_re, function(i) paste0(" + XS_RE_", i, "[n] * beta_S_RE_", i, "[s-1]'")), collapse = '')),
      ";"
    ),
    "        }"
  )
  b_else <- F
  ZD_comb <- unique(PSobject$SZDG_table[, c("Z", "D")])
  for (r in 1:nrow(ZD_comb)) {
    z <- ZD_comb[r, "Z"]
    d <- ZD_comb[r, "D"]
    SG_table <- PSobject$SZDG_table[PSobject$SZDG_table[, "Z"] == z & PSobject$SZDG_table[, "D"] == d, , drop = F]
    model_output_lines <- c(
      model_output_lines,
      paste0("        ", ifelse(b_else, "else ", ""), "if (Z[n] == ", z, " && D[n] == ", d, ")"),
      paste0("            length = ", nrow(SG_table), ";")
    )
    b_else <- T
  }
  model_output_lines <- c(
    model_output_lines,
    "        {",
    "            real log_l[length];"
  )
  b_else <- F
  for (r in 1:nrow(ZD_comb)) {
    z <- ZD_comb[r, "Z"]
    d <- ZD_comb[r, "D"]
    SG_table <- PSobject$SZDG_table[PSobject$SZDG_table[, "Z"] == z & PSobject$SZDG_table[, "D"] == d, , drop = F]
    model_output_lines <- c(
      model_output_lines,
      paste0("            ", ifelse(b_else, "else ", ""), "if (Z[n] == ", z, " && D[n] == ", d, ") {"),
      paste0("                // Z:", z, " D:", d, " S:", paste(SG_table[, "S"], collapse = "/"))
    )
    number <- 1
    for (sg in 1:nrow(SG_table)) {
      s <- SG_table[sg, 'S']
      g <- SG_table[sg, 'G']
      str_param <- paste0(
        ", ",
        paste(sapply(parameter_list, function(param) paste0(param$name, "[", g, "]")), collapse = ", "),
        ifelse (Y_type == "survival", ", delta[n]", "")
      )
      model_output_lines <- c(
        model_output_lines,
        paste0("                log_l[", number, "] = log_prob[", s + 1, "] + ", func_family, "(Y[n] | ",
               func_link, "(", "XG[n] * beta_G[", g, "]'", 
               ifelse(Y_re == 0, "", paste(
                 1:Y_re, function(i) paste0(" + XG_RE_", i, "[n] * beta_G_RE_", i, "[", g, "]'")
               )),
               ")", str_param, ");"
        )
      )
      number <- number + 1
    }
    model_output_lines <- c(
      model_output_lines,
      "            }"
    )
    b_else <- T
  }
  model_output_lines <- c(
    model_output_lines,
    "            target += log_sum_exp(log_l) - log_sum_exp(log_prob);",
    "        }",
    "    }",
    "}"
  )
  
  ## generated_quantities
  generated_quantities_output_lines <- c(
    "generated quantities {",
    paste0("    vector[", SZDG_max['S'] + 1, "] strata_prob; // the probability of being in each stratum"),
    ifelse(Y_type == "survival",
           paste0("    matrix[", SZDG_max['G'], ", T] mean_surv_prob; // mean survival probability"),
           paste0("    vector[", SZDG_max['G'], "] mean_effect; // mean response")),
    "    {",
    paste0("        matrix[N, ", SZDG_max['S'] + 1, "] log_prob;")
  )
  if(Y_type == "survival"){
    generated_quantities_output_lines <- c(
      generated_quantities_output_lines,
      paste0("        matrix[", SZDG_max['G'], ", T] numer;"),
      paste0("        matrix[N, T] expected_surv_prob[", SZDG_max['G'], "];"),
      "        for (i in 1:N)",
      paste0("            for (j in 1:", SZDG_max['G'], ")"),
      "                for (t in 1:T) {",
      paste0("                    real mu = ", func_link, "(XG[i] * beta_G[j]'",
             ifelse(Y_re == 0, "", 
                    paste(sapply(1:Y_re, function(i) paste0(" + XG_RE_", i, "[i] * beta_G_RE_", i, "[j]'")))),
             ");"
      ),
      ifelse(Y_family == "survival_Cox", 
             "                    expected_surv_prob[j][i,t] = exp(-exp(mu) * pow(time[t], exp(theta[j])));", ""),
      "                }"
    )
  } else {
    generated_quantities_output_lines <- c(
      generated_quantities_output_lines,
      paste0("        vector[", SZDG_max['G'], "] numer;"),
      paste0('        matrix[N, ', SZDG_max['G'], "] expected_mean;"),
      "        for (i in 1:N)",
      paste0("            for (j in 1:", SZDG_max['G'], ")"),
      paste0("                expected_mean[i, j] = ", func_link, "(XG[i] * beta_G[j]'", 
             ifelse(Y_re == 0, "", 
                    paste(sapply(1:Y_re, function(i) paste0(" + XG_RE_", i, "[i] * beta_G_RE_", i, "[j]'")))),
             ");"
      )
    )
  }
  generated_quantities_output_lines <- c(
    generated_quantities_output_lines,
    "        log_prob[:, 1] = rep_vector(0, N);",
    paste0("        log_prob[:, 2:", SZDG_max['S'] + 1, "] = XS * beta_S'",
           ifelse(S_re == 0, "", 
                  paste(sapply(1:S_re, function(i) paste0(" + XS_RE_", i, " * beta_S_RE_", i, "'")))),
           ";"
    ),
    "        for (n in 1:N) {",
    "            log_prob[n] -= log_sum_exp(log_prob[n]);",
    "        }",
    paste0("        for (s in 1:", SZDG_max['S'] + 1, ") strata_prob[s] = mean(exp(log_prob[:, s]));")
  )
  if (Y_type == "survival") {
    generated_quantities_output_lines <- c(
      generated_quantities_output_lines,
      paste0("        for (g in 1:", SZDG_max['G'], ") {"),
      "            for (t in 1:T) {",
      "                numer[g, t] = mean(expected_surv_prob[g][:, t] .* exp(log_prob[:, S[g]]));",
      "                mean_surv_prob[g, t] = numer[g, t] / strata_prob[S[g]];",
      "            }",
      "        }"
    )
  } else {
    generated_quantities_output_lines <- c(
      generated_quantities_output_lines,
      paste0("        for (g in 1:", SZDG_max['G'], ") {"),
      "            numer[g] = mean(expected_mean[:, g] .* exp(log_prob[:, S[g]]));",
      "            mean_effect[g] = numer[g] / strata_prob[S[g]];",
      "        }"
    )
  }
  generated_quantities_output_lines <- c(
    generated_quantities_output_lines,
    "    }",
    "}"
  )
  
  comment_lines <- c(
    "// Automatically generated by PStrata 0.0.1",
    paste0("// at ", Sys.time(), " ", Sys.timezone())
  )
  
  full_output <- c(
    comment_lines, " ",
    func_output_lines, " ",
    stan_data_output_lines, " ",
    transformed_data_output_lines, " ",
    parameters_output_lines, " ",
    transformed_parameters_output_lines, " ",
    model_output_lines, " ",
    generated_quantities_output_lines
  )
  
  full_text <- paste0(full_output, collapse = "\n")
  if (!is.null(filename)) {
    if (!stringr::str_ends(filename, ".stan"))
      filename <- paste0(filename, ".stan")
    cat(full_text, file = filename)
  }
  return (full_text)
}
