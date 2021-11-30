generate_names <- function(class, group, name, dim = NULL) {
  tmp_str <- paste0(c(class, group, name), collapse = '_')
  if (!is.null(dim)){
    if (dim > 0)
      return (paste(tmp_str, "[", 1:dim, "]", sep = ''))
    else
      return (c())
  }
  else
    return (tmp_str)
}

PSSampling <- function(PSobject, model_name = "unnamed", ...) {
  write.pso(PSobject, paste0(model_name, ".pso"))
  to_stan(model_name)
  stanfit <- rstan::stan(paste0(model_name, ".stan"), 
                         data = PSobject$data, 
                         ...
  )
  post_samples <- do.call(cbind, rstan::extract(stanfit))
  param_names <- c()
  for (param in PSobject$parameter_list) {
    param_names <- c(param_names, 
                     generate_names(param$class, param$group, param$name, param$dim))
  }
  param_names <- c(param_names, '__lp')
  colnames(post_samples) <- param_names
  return (structure(
    list(
      stanfit = stanfit,
      post_samples = post_samples
    ),
    class = "PSsample"
  ))
}

PSSampleEx <- function(PSobject, PSsample) {
  exponential <- function(log_prob_array) {
    max_log_prob <- apply(log_prob_array, c(2, 3), max)
    return (exp(sweep(log_prob_array, c(2, 3), max_log_prob)))
  }
  standardize <- function(prob_array) {
    sum_prob <- apply(prob_array, c(2, 3), sum)
    return (sweep(prob_array, c(2, 3), sum_prob, FUN = '/'))
  }
  
  dimension <- c(length(PSobject$PSsettings$strata), 
                 PSobject$PScovs$ncount, nrow(PSsample$post_samples))
  log_post_prob_unadjusted <- array(
    0, dim = dimension, dimnames = list(PSobject$PSsettings$strata, NULL, NULL)
  ) # log_probability
  for (stratum in PSobject$PSsettings$strata[-1]) {
    value <- 0
    if (generate_names('S', stratum, 'Intrcpt') %in% colnames(PSsample$post_samples))
      intercept <- PSsample$post_samples[, 
                            generate_names('S', stratum, 'Intrcpt'),
                            drop = F]
    else
      intercept <- 0
    value <- value + rep(1, each = PSobject$PScovs$ncount) %*% t(intercept)
    if (PSobject$PScovs$S.dim > 0){
      pars <- PSsample$post_samples[, generate_names(
        'S', stratum, 'Coef', PSobject$PScovs$S.dim
      )]
      value <- value + PSobject$data$XS %*% t(pars)
    }
    log_post_prob_unadjusted[stratum, , ] <- value
  }
  
  # consistency matrix
  cs_mat <- array(
    1, dim = dimension, dimnames = list(PSobject$PSsettings$strata, NULL, NULL)
  )
  for (stratum in PSobject$PSsettings$strata) {
    cs_mat[stratum, , ] <- ifelse(
      substring(stratum, PSobject$data$Z+1, PSobject$data$Z+1) == PSobject$data$D, 1, 0
    )
  }
  
  post_prob_unadjusted <- standardize(exponential(log_post_prob_unadjusted))
  post_prob_consistent <- standardize(post_prob_unadjusted * cs_mat)
  post_stratum_draw <- apply(
    post_prob_consistent, c(2,3), 
    function(p) 
      PSobject$PSsettings$strata[(1:dim(post_prob_consistent)[1]) %*% rmultinom(1,1,p)]
  )
  
  post_outcome_mean <- array(
    0, dim = c(dimension[1], 2, dimension[-1]),
    dimnames = list(PSobject$PSsettings$strata, c("0", "1"), NULL, NULL)
  )
  theta <- array(
    0, dim = c(dimension[1], 2, dimension[3]),
    dimnames = list(PSobject$PSsettings$strata, c("0", "1"), NULL)
  )
  for (stratum in PSobject$PSsettings$strata) {
    for (z in c("0", "1")){
      zz <- ifelse(stratum %in% PSobject$PSsettings$ER, "0", z)
      value <- 0
      if (generate_names('Y', c(stratum, zz), 'Intrcpt') %in% colnames(PSsample$post_samples))
        intercept <- PSsample$post_samples[, 
                                           generate_names('Y', c(stratum, zz), 'Intrcpt'), 
                                           drop = F]
      else
        intercept <- 0
      value <- value + rep(1, each = PSobject$PScovs$ncount) %*% t(intercept)
      if (PSobject$PScovs$Y.dim > 0){
        pars <- PSsample$post_samples[, generate_names('Y', c(stratum, zz), 'Coef', PSobject$PScovs$Y.dim)]
        value <- value + PSobject$data$XY %*% t(pars)
      }
      if (PSobject$PSsettings$Y.family$family == "survival") {
        theta[stratum, z, ] <- PSsample$post_samples[, generate_names('Y', c(stratum, zz), 'Theta')]
      }
      post_outcome_mean[stratum, z, , ] <- value
    }
  }
  if (PSobject$PSsettings$Y.family$family != "survival") {
    post_outcome <- list(
      outcome = PSobject$PSsettings$Y.family$linkinv(post_outcome_mean)
    )
  }
  else {
    post_outcome <- list(
      outcome = post_outcome_mean,
      theta = theta
    )
  }

  return (structure(
    list(
      PSobject = PSobject,
      PSsample = PSsample,
      prob_unadj = post_prob_unadjusted,
      prob_const = post_prob_consistent,
      consistency = cs_mat,
      stratum_draw = post_stratum_draw,
      outcome = post_outcome
    ),
    class = "PSSampleEx",
    survival = PSobject$PSsettings$Y.family$family == "survival"
  ))
}

PSSummary <- function(PSsampleEx){
  if (attr(PSsampleEx, "survival"))
    return (PSSummary.survival(PSsampleEx))
  else
    return (PSSummary.nonsurvival(PSsampleEx))
}

PSSummary.nonsurvival <- function(PSsampleEx) {
  prob <- broadcast(extend(PSsampleEx$prob_const, 2), dim(PSsampleEx$outcome$outcome))
  outcome <- apply(PSsampleEx$outcome$outcome * prob, c(1, 2, 4), mean) / apply(prob, c(1, 2, 4), mean)
  causal_effect <- outcome[, "1", ] - outcome[, "0", ]

  return (structure(list(
    PSsampleEx = PSsampleEx,
    Overall_0 = outcome[, "0", ],
    Overall_1 = outcome[, "1", ],
    Causal_Effect = causal_effect
  ), class = "PSSummary.nonsurvival"))
}

PSSummary.survival <- function(PSsampleEx) {
  hazard_at <- function(time) {
    full_dim <- dim(PSsampleEx$outcome$outcome)
    theta <- broadcast(extend(PSsampleEx$outcome$theta, 3), full_dim)
    mu <- PSsampleEx$outcome$outcome
    log_t <- log(time)
    log_S <- -exp(mu) * exp(log_t * exp(theta))
    log_lambda <- theta + mu + (exp(theta) - 1) * log_t
    prob <- broadcast(extend(PSsampleEx$prob_const, c(2)), full_dim)
    survival_curve <- apply(exp(log_S) * prob, c(1,2,4), mean) / apply(prob, c(1,2,4), mean)
    max_log <- max(log_S, log_S + log_lambda)
    hazard <- apply(exp(log_S + log_lambda - max_log) * prob, c(1,2,4), mean) / 
      apply(exp(log_S - max_log) * prob, c(1,2,4), mean)
    hazard_ratio <- hazard[, "1",] / hazard[, "0",]
    return (list(
      survival_curve = survival_curve,
      hazard_curve = hazard,
      hazard_ratio = hazard_ratio
    ))
  }

  return (structure(list(
    PSsampleEx = PSsampleEx,
    hazard_at = hazard_at
  ), class = "PSSummary.survival"))
}

print.PSSummary.nonsurvival <- function(PSsummary) {
  cat("Estimated Proportion from Each Stratum:\n")
  prob_samples <- apply(PSsummary$PSsampleEx$prob_const, c(1, 3), mean)
  prob_summary <- apply(prob_samples, 1, function(x) {
    c(mean = mean(x), sd = sd(x), 
      lwr = quantile(x, 0.025), upr = quantile(x, 0.975))
  })
  print(prob_summary)
  cat('\n')
  cat("Causal Effect of Principal Strata:\n")
  outcome_summary <- apply(PSsummary$Causal_Effect, 1, function(x) {
    c(mean = mean(x), sd = sd(x), 
      lwr = quantile(x, 0.025), upr = quantile(x, 0.975))
  })
  print(outcome_summary)
  cat('\n')
}

print.PSSummary.survival <- function(PSsummary) {
  cat("Estimated Proportion from Each Stratum:\n")
  prob_samples <- apply(PSsummary$PSsampleEx$prob_const, c(1, 3), mean)
  prob_summary <- apply(prob_samples, 1, function(x) {
    c(mean = mean(x), sd = sd(x), 
      lwr = quantile(x, 0.025), upr = quantile(x, 0.975))
  })
  print(prob_summary)
  cat('\n')
  cat("Causal Effect of Principal Strata:\n")
  cat("Call summary() with time parameter to get hazard ratio.")
  cat('\n')
}

summary.PSSummary.nonsurvival <- function(PSsummary){
  func_summarize <- function(x) {
    c(mean = mean(x), 
      se_mean = ifelse(abs(sd(x)) < .Machine$double.eps, 
                       0, 
                       sd(x) / unname(sqrt(coda::effectiveSize(x)))), 
      sd = sd(x), 
      quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))
  }
  parameter <- rstan::summary(PSsummary$PSsampleEx$PSsample$stanfit)$summary
  
  # proportion
  prob_samples <- apply(PSsummary$PSsampleEx$prob_const, c(1, 3), mean)
  prob_summary <- apply(prob_samples, 1, func_summarize)
  proportion <- t(prob_summary)
  
  # outcome
  out_summary_0 <- apply(PSsummary$Overall_0, 1, func_summarize)
  out_summary_1 <- apply(PSsummary$Overall_1, 1, func_summarize)
  effect <- apply(PSsummary$Causal_Effect, 1, func_summarize)
  outcome_0 <- t(out_summary_0)
  outcome_1 <- t(out_summary_1)
  effect <- t(effect)
  
  return (structure(list(
    parameter = parameter,
    proportion = proportion,
    outcome = list(outcome_0 = outcome_0, outcome_1 = outcome_1),
    effect = effect
  ), class = "summary.PSSummary.nonsurvival"))
}

summary.PSSummary.survival <- function(PSsummary, time = 1){
  func_summarize <- function(x) {
    c(mean = mean(x), 
      se_mean = ifelse(abs(sd(x)) < .Machine$double.eps, 
                       0, 
                       sd(x) / unname(sqrt(coda::effectiveSize(x)))), 
      sd = sd(x), 
      quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))
  }
  parameter <- rstan::summary(PSsummary$PSsampleEx$PSsample$stanfit)$summary
  
  # proportion
  prob_samples <- apply(PSsummary$PSsampleEx$prob_const, c(1, 3), mean)
  prob_summary <- apply(prob_samples, 1, func_summarize)
  proportion <- t(prob_summary)
  
  # outcome
  strata <- PSsummary$PSsampleEx$PSobject$PSsettings$strata
  survival_result <- pbapply::pblapply(time, PSsummary$hazard_at)
  dim_full <- c(dim(survival_result[[1]]$survival_curve), length(time))
  dimnames_full <- c(dimnames(survival_result[[1]]$survival_curve), list(time))
  survival_curve_full <- array(NA, dim = dim_full, dimnames = dimnames_full)
  hazard_curve_full <- array(NA, dim = dim_full, dimnames = dimnames_full)
  hazard_ratio_full <- array(NA, dim = dim_full[-2], dimnames = dimnames_full[-2])
  for (i in 1:length(time)) {
    survival_curve_full[,,,i] <- survival_result[[i]]$survival_curve
    hazard_curve_full[,,,i] <- survival_result[[i]]$hazard_curve
    hazard_ratio_full[,,i] <- survival_result[[i]]$hazard_ratio
  }
  
  
  hazard_ratio <- hazard_curve <- survival_curve <- list()
  
  dimension <- dim(survival_curve_full)[3:4]
  for (stratum in strata){
    for (trt in c("0", "1")) {
      surv_func <- t(apply(
        array(survival_curve_full[stratum, trt,,], dim = dimension),
        2, func_summarize
      ))
      rownames(surv_func) <- time
      
      hzrd_func <- t(apply(
        array(hazard_curve_full[stratum, trt,,], dim = dimension),
        2, func_summarize
      ))
      rownames(hzrd_func) <- time
      
      survival_curve[[stratum]][[trt]] <- surv_func
      hazard_curve[[stratum]][[trt]] <- hzrd_func
    }
    
    hzrd_ratio <- t(apply(
      array(hazard_ratio_full[stratum,,], dim = dimension),
      2, func_summarize
    ))
    rownames(hzrd_ratio) <- time
    hazard_ratio[[stratum]] <- hzrd_ratio
  }
  
  return (structure(list(
    time_points = time,
    parameter = parameter,
    proportion = proportion,
    survival_curve = survival_curve,
    hazard_curve = hazard_curve,
    hazard_ratio = hazard_ratio
  ), class = "summary.PSSummary.survival"))
}

print.summary.PSSummary.nonsurvival <- function(summary.PSsummary) {
  
  cat("Posterior estimate of parameters:\n")
  print(summary.PSsummary$parameter)
  cat("\n")
  
  cat("Estimated Proportion from Each Stratum:\n")
  print(summary.PSsummary$proportion)
  cat("\n")
  
  cat("Estimated Average Outcome:\n")
  cat("-- Treatment (Z = 1) --\n")
  print(summary.PSsummary$outcome$outcome_1)
  cat("\n")
  cat("-- Control (Z = 0) --\n")
  print(summary.PSsummary$outcome$outcome_0)
  cat("\n")
  
  cat("Estimated Treatment Effect (Control as Reference):\n")
  print(summary.PSsummary$effect)
  cat("\n")
}

print.summary.PSSummary.survival <- function(summary.PSsummary) {
  func_summarize <- function(x) {
    c(mean = mean(x, na.rm = T), 
      se_mean = ifelse(abs(sd(x, na.rm = T)) < .Machine$double.eps, 
                       0, 
                       sd(x, na.rm = T) / unname(sqrt(coda::effectiveSize(x)))), 
      sd = sd(x, na.rm = T), 
      quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = T))
  }
  cat("Posterior estimate of parameters:\n")
  print(summary.PSsummary$parameter)
  cat("\n")
  
  cat("Estimated Proportion from Each Stratum:\n")
  print(summary.PSsummary$proportion)
  cat("\n")
  
  cat("Estimated Survival Function (Mean):\n")
  for (stratum in rownames(summary.PSsummary$proportion)) {
    cat("Stratum ", stratum, ":\n", sep = "")
    cat("-- Treatment (Z = 1) --\n")
    print(summary.PSsummary$survival_curve[[stratum]][["1"]])
    cat("\n")
    cat("-- Control (Z = 0) --\n")
    print(summary.PSsummary$survival_curve[[stratum]][["0"]])
    cat("\n")
  }
  
  cat("Estimated Hazard Function (Mean):\n")
  for (stratum in rownames(summary.PSsummary$proportion)) {
    cat("Stratum ", stratum, ":\n", sep = "")
    cat("-- Treatment (Z = 1) --\n")
    print(summary.PSsummary$hazard_curve[[stratum]][["1"]])
    cat("\n")
    cat("-- Control (Z = 0) --\n")
    print(summary.PSsummary$hazard_curve[[stratum]][["0"]])
    cat("\n")
  }
  
  cat("Estimated Hazard Ratio (Control as Reference):\n")
  for (stratum in rownames(summary.PSsummary$proportion)) {
    cat("Stratum ", stratum, ":\n", sep = "")
    print(summary.PSsummary$hazard_ratio[[stratum]])
    cat("\n")
  }
}

PStrata <- function(S.formula, Y.formula, Y.family, data, monotonicity = "default", ER = c(), trunc = FALSE, 
               prior_intercept = prior_uniform(),
               prior_coefficient = prior_normal(),
               prior_sigma = prior_inv_gamma(),
               prior_alpha = prior_inv_gamma(),
               prior_lambda = prior_inv_gamma(),
               prior_theta = prior_normal(),
               ...){
  PSobject <- PSObject(
    S.formula, Y.formula, Y.family, data, monotonicity, ER, trunc,
    prior_intercept, prior_coefficient, prior_sigma,
    prior_alpha, prior_lambda, prior_theta
  )
  PSsample <- PSSampling(PSobject, ...)
  PSsampleEx <- PSSampleEx(PSobject, PSsample)
  PSsummary <- PSSummary(PSsampleEx)
  res <- list(
    PSobject = PSobject,
    PSsummary = PSsummary
  )
  class(res) <- "PStrata"
  return (res)
}

plot.PSSummary.nonsurvival <- function(PSsummary) {
  get_long_df <- function(table, var) {
    df <- tidyr::pivot_longer(as.data.frame(t(table)), everything())
    df$var <- var
    return (df)
  }
  probability <- apply(PSsummary$PSsampleEx$prob_const, c(1, 3), mean)
  long_prob <- get_long_df(probability, var = "Probability")
  long_O1 <- get_long_df(PSsummary$Overall_1, var = "Overall_1")
  long_O0 <- get_long_df(PSsummary$Overall_0, var = "Overall_0")
  long_CE <- get_long_df(PSsummary$Causal_Effect, var = "Causal Effect")
  table_1 <- dplyr::bind_rows(long_O1, long_O0)
  p0 <- ggplot(long_prob) + 
    geom_density(aes(x = value, color = name)) +
    geom_histogram(aes(x = value, y = after_stat(density), fill = name), alpha = 0.3) +
    guides(color = "none", fill = "none") +
    ggtitle("Probability") + xlab("probability")
  p1 <- ggplot(table_1) + 
    geom_density(aes(x = value, color = var)) +
    geom_histogram(aes(x = value, y = after_stat(density), fill = var), alpha = 0.3) +
    facet_wrap(~name) + ggtitle("Mean Effect")
  p2 <- ggplot(long_CE) + 
    geom_histogram(aes(x = value, fill = name), alpha = 0.8) +
    ggtitle("Average Treatment Effect")
  (p0 + p2) / p1
}

plot.summary.PSSummary.survival <- function(summary.PSsummary, PSsummary) {
  get_long_df <- function(table, var) {
    df <- tidyr::pivot_longer(as.data.frame(t(table)), everything())
    df$var <- var
    return (df)
  }
  probability <- apply(PSsummary$PSsampleEx$prob_const, c(1, 3), mean)
  survival_result <- summary.PSsummary
  long_prob <- get_long_df(probability, var = "Probability")
  strata <- PSsummary$PSsampleEx$PSobject$PSsettings$strata
  
  time <- summary.PSsummary$time_points
  table_list_hazard_ratio <- list()
  for (stratum in strata){
    tmp <- survival_result$hazard_ratio[[stratum]]
    table_list_hazard_ratio <- c(
      table_list_hazard_ratio,
      list(data.frame(
        stratum = stratum,
        time = time, 
        est = tmp[, "50%"],
        lwr = tmp[, "2.5%"],
        upr = tmp[, "97.5%"]
      )
      ))
  }
  long_table_hazard_ratio <- do.call(dplyr::bind_rows, table_list_hazard_ratio)
  
  table_list_hazard_curve <- list()
  for (stratum in strata){
    for (trt in c("0", "1")){
      tmp <- survival_result$hazard_curve[[stratum]][[trt]]
      table_list_hazard_curve <- c(
        table_list_hazard_curve,
        list(data.frame(
          stratum = stratum,
          treatment = trt,
          time = time, 
          est = tmp[, "50%"],
          lwr = tmp[, "2.5%"],
          upr = tmp[, "97.5%"]
        )
        ))
    }
  }
  long_table_hazard_curve <- do.call(dplyr::bind_rows, table_list_hazard_curve)
  
  table_list_survival_curve <- list()
  for (stratum in strata){
    for (trt in c("0", "1")){
      tmp <- survival_result$survival_curve[[stratum]][[trt]]
      table_list_survival_curve <- c(
        table_list_survival_curve,
        list(data.frame(
          stratum = stratum,
          treatment = trt,
          time = time, 
          est = tmp[, "50%"],
          lwr = tmp[, "2.5%"],
          upr = tmp[, "97.5%"]
        )
        ))
    }
  }
  long_table_survival_curve <- do.call(dplyr::bind_rows, table_list_survival_curve)
  
  p0 <- ggplot(long_prob) + 
    geom_density(aes(x = value, color = name)) +
    geom_histogram(aes(x = value, y = after_stat(density), fill = name), alpha = 0.3) +
    guides(color = "none", fill = "none") +
    ggtitle("Probability") + xlab("probability")
  p1 <- ggplot(long_table_hazard_ratio) + 
    geom_line(aes(x = time, y = est, group = stratum, color = stratum)) +
    geom_ribbon(aes(x = time, ymin = lwr, ymax = upr, fill = stratum), alpha = 0.3) +
    scale_y_log10() +
    ggtitle("Hazard Ratio")
  p2 <- ggplot(long_table_hazard_curve) + 
    geom_line(aes(x = time, y = est, group = treatment, color = treatment)) +
    geom_ribbon(aes(x = time, ymin = lwr, ymax = upr, fill = treatment), alpha = 0.3) +
    facet_wrap(~stratum, scale = "free_y") +
    ggtitle("Hazard Curve")
  p3 <- ggplot(long_table_survival_curve) + 
    geom_line(aes(x = time, y = est, group = treatment, color = treatment)) +
    geom_ribbon(aes(x = time, ymin = lwr, ymax = upr, fill = treatment), alpha = 0.3) +
    facet_wrap(~stratum, scale = "free_y") +
    ggtitle("Survival Curve")
  (p0 + p1) / p2 / p3
}

plot.PSSummary.survival <- function(PSsummary, time = 1) {
  get_long_df <- function(table, var) {
    df <- tidyr::pivot_longer(as.data.frame(t(table)), everything())
    df$var <- var
    return (df)
  }
  probability <- apply(PSsummary$PSsampleEx$prob_const, c(1, 3), mean)
  survival_result <- summary.PSSummary.survival(PSsummary, time)
  long_prob <- get_long_df(probability, var = "Probability")
  strata <- PSsummary$PSsampleEx$PSobject$PSsettings$strata
  
  table_list_hazard_ratio <- list()
  for (stratum in strata){
    tmp <- survival_result$hazard_ratio[[stratum]]
    table_list_hazard_ratio <- c(
      table_list_hazard_ratio,
      list(data.frame(
        stratum = stratum,
        time = time, 
        est = tmp[, "50%"],
        lwr = tmp[, "2.5%"],
        upr = tmp[, "97.5%"]
      )
    ))
  }
  long_table_hazard_ratio <- do.call(dplyr::bind_rows, table_list_hazard_ratio)
  
  table_list_hazard_curve <- list()
  for (stratum in strata){
    for (trt in c("0", "1")){
      tmp <- survival_result$hazard_curve[[stratum]][[trt]]
      table_list_hazard_curve <- c(
        table_list_hazard_curve,
        list(data.frame(
          stratum = stratum,
          treatment = trt,
          time = time, 
          est = tmp[, "50%"],
          lwr = tmp[, "2.5%"],
          upr = tmp[, "97.5%"]
        )
        ))
    }
  }
  long_table_hazard_curve <- do.call(dplyr::bind_rows, table_list_hazard_curve)
  
  table_list_survival_curve <- list()
  for (stratum in strata){
    for (trt in c("0", "1")){
      tmp <- survival_result$survival_curve[[stratum]][[trt]]
      table_list_survival_curve <- c(
        table_list_survival_curve,
        list(data.frame(
          stratum = stratum,
          treatment = trt,
          time = time, 
          est = tmp[, "50%"],
          lwr = tmp[, "2.5%"],
          upr = tmp[, "97.5%"]
        )
        ))
    }
  }
  long_table_survival_curve <- do.call(dplyr::bind_rows, table_list_survival_curve)
  
  p0 <- ggplot(long_prob) + 
    geom_density(aes(x = value, color = name)) +
    geom_histogram(aes(x = value, y = after_stat(density), fill = name), alpha = 0.3) +
    guides(color = "none", fill = "none") +
    ggtitle("Probability") + xlab("probability")
  p1 <- ggplot(long_table_hazard_ratio) + 
    geom_line(aes(x = time, y = est, group = stratum, color = stratum)) +
    geom_ribbon(aes(x = time, ymin = lwr, ymax = upr, fill = stratum), alpha = 0.3) +
    ggtitle("Hazard Ratio")
  p2 <- ggplot(long_table_hazard_curve) + 
    geom_line(aes(x = time, y = est, group = treatment, color = treatment)) +
    geom_ribbon(aes(x = time, ymin = lwr, ymax = upr, fill = treatment), alpha = 0.3) +
    facet_wrap(~stratum) +
    ggtitle("Hazard Curve")
  p3 <- ggplot(long_table_survival_curve) + 
    geom_line(aes(x = time, y = est, group = treatment, color = treatment)) +
    geom_ribbon(aes(x = time, ymin = lwr, ymax = upr, fill = treatment), alpha = 0.3) +
    facet_wrap(~stratum) +
    ggtitle("Survival Curve")
  (p0 + p1) / p2 / p3
}

print.PStrata <- function(pstrata) {
  print(pstrata$PSsummary)
}

summary.PStrata <- function(pstrata) {
  summary(pstrata$PSsummary)
}


plot_prob_prob_one <- function(num_draw, post_prob, post_stratum_samples) {
  post_stratum <- post_stratum_samples[, num_draw]
  prob <- post_prob[,, num_draw]
  data <- bind_cols(as.data.frame(t(prob)), stratum = post_stratum)
  colnames(data) <- c(obj$strata, "stratum")
  data$stratum <- obj$strata[data$stratum]
  data_plot <- pivot_longer(data, -stratum)
  colnames(data_plot) <- c("sample", "stratum", "prob")
  ggplot(data_plot) + 
    geom_histogram(aes(prob, y = ..density.., fill = stratum), alpha = 0.4) + 
    geom_density(aes(prob, color = stratum)) + 
    facet_wrap(~sample)
}

plot.PStrata <- function(res, ...){
  plot(res$outcome, ...)
}
