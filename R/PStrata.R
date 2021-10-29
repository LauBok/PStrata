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
        theta[stratum, z, ] <- PSobject$post_samples[, generate_names('Y', c(stratum, zz), 'Theta')]
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
  Y0 <- PSsampleEx$outcome$outcome[, "0", ,]
  Y1 <- PSsampleEx$outcome$outcome[, "1", ,]
  sum_Y0 <- apply(Y0 * PSsampleEx$prob_const, c(1, 3), sum)
  sum_Y1 <- apply(Y1 * PSsampleEx$prob_const, c(1, 3), sum)
  sum_post_prob <- apply(PSsampleEx$prob_const, c(1, 3), sum)
  mean_Y0 <- sum_Y0 / sum_post_prob
  mean_Y1 <- sum_Y1 / sum_post_prob
  causal_effect <- mean_Y1 - mean_Y0
  return (structure(list(
    Overall_0 = mean_Y0,
    Overall_1 = mean_Y1,
    Causal_Effect = causal_effect
  ), class = "PSsummary.nonsurvival"))
}

PSSummary.survival <- function(PSsampleEx) {
  Y0 <- PSsampleEx$outcome$outcome[, "0", ,]
  Y1 <- PSsampleEx$outcome$outcome[, "1", ,]
  theta0 <- PSsampleEx$outcome$theta[, "0", ]
  theta1 <- PSsampleEx$outcome$theta[, "1", ]
  hazard_at <- function (time) {
    # log(wi)
    log_w0 <- sweep(-exp(Y0), c(1, 3), time^(exp(theta0)), FUN = '*')
    log_w1 <- sweep(-exp(Y1), c(1, 3), time^(exp(theta1)), FUN = '*')
    log_w0 <- sweep(log_w0, c(1, 3), apply(log_w0, c(1, 3), max))
    log_w1 <- sweep(log_w1, c(1, 3), apply(log_w1, c(1, 3), max))
    wp0 <- exp(log_w0) * PSsampleEx$prob_const
    wp1 <- exp(log_w1) * PSsampleEx$prob_const
    
    hzrd0_base <- exp(sweep(Y0, c(1,3), theta0, FUN = '+'))
    hzrd1_base <- exp(sweep(Y1, c(1,3), theta1, FUN = '+'))
    hzrd0 <- sweep(hzrd0_base, c(1, 3), time^(exp(theta0) - 1), FUN = '*')
    hzrd1 <- sweep(hzrd1_base, c(1, 3), time^(exp(theta1) - 1), FUN = '*')
    sum_hzrd0 <- apply(hzrd0 * wp0, c(1, 3), sum)
    sum_hzrd1 <- apply(hzrd1 * wp1, c(1, 3), sum)
    sum_wp0 <- apply(wp0, c(1, 3), sum)
    sum_wp1 <- apply(wp1, c(1, 3), sum)
    mean_hzrd0 <- sum_hzrd0 / sum_wp0
    mean_hzrd1 <- sum_hzrd1 / sum_wp1
    hazard_ratio <- mean_hzrd0 / mean_hzrd1
    return (list(
      Hazard_0 = mean_hzrd0,
      Hazard_1 = mean_hzrd1,
      Hazard_Ratio = hazard_ratio
    ))
  }
  return (structure(list(
    hazard_at = hazard_at
  ), class = "PSSummary.survival"))
}

PStrata <- function(S.formula, Y.formula, Y.family, data, monotonicity = "default", ER = c(), trunc = FALSE, 
               prior_intercept = prior_uniform(),
               prior_coefficient = prior_normal(),
               prior_sigma = prior_inv_gamma(),
               prior_alpha = prior_inv_gamma(),
               prior_lambda = prior_inv_gamma(),
               prior_theta = prior_normal(),
               ...){
  PSobj <- PSObject(
    S.formula, Y.formula, Y.family, data, monotonicity, ER,
    prior_intercept, prior_coefficient, prior_sigma,
    prior_alpha, prior_lambda, prior_theta
  )
  write.pso(PSobj, paste0(model_name, ".pso"))
  to_stan(model_name)
  df <- get_data(obj, data)
  stanfit <- rstan::stan(paste0(model_name, ".stan"), data = df, 
                         ...
  )
  post_samples <- get_posterior_samples(stanfit, obj)
  post_probability_raw <- post_prob_raw(obj, post_samples, df)
  consistency_mat <- consistency(obj, post_samples, df)
  post_probability <- post_prob(post_probability_raw, consistency_mat)
  post_outcome <- post_Y(obj, post_samples, df)
  if (obj$Y.family$family == "survival") {
    outcome <- summarize_result_survival(post_probability_raw, post_probability, post_outcome)
  }
  else
    outcome <- summarize_result(post_probability_raw, post_probability, post_outcome)
  
  res <- list(
    data = df,
    stanfit = stanfit,
    PSobject = obj,
    pso_code = pso_code,
    stan_code = stan_code,
    post_samples = post_samples,
    outcome = outcome
  )
  class(res) <- "PStrata"
  return (res)
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

summarize_result <- function(post_prob_raw, post_prob, post_Y) {
  Y0 <- post_Y[, "0", ,]
  Y1 <- post_Y[, "1", ,]
  sum_Y0 <- apply(Y0 * post_prob, c(1, 3), sum)
  sum_Y1 <- apply(Y1 * post_prob, c(1, 3), sum)
  sum_post_prob <- apply(post_prob, c(1, 3), sum)
  mean_Y0 <- sum_Y0 / sum_post_prob
  mean_Y1 <- sum_Y1 / sum_post_prob
  causal_effect <- mean_Y1 - mean_Y0
  return (structure(list(
    Prob = post_prob_raw,
    Overall_0 = mean_Y0,
    Overall_1 = mean_Y1,
    Causal_Effect = causal_effect
  ), class = "post_outcome"))
}

summarize_result_survival <- function(post_prob_raw, post_prob, post_Y) {
  Y0 <- post_Y$mean[, "0", ,]
  Y1 <- post_Y$mean[, "1", ,]
  theta0 <- post_Y$theta[, "0", ]
  theta1 <- post_Y$theta[, "1", ]
  f <- function (t) {
    # log(wi)
    log_w0 <- sweep(-exp(Y0), c(1, 3), t^(exp(theta0)), FUN = '*')
    log_w1 <- sweep(-exp(Y1), c(1, 3), t^(exp(theta1)), FUN = '*')
    log_w0 <- sweep(log_w0, c(1, 3), apply(log_w0, c(1, 3), max))
    log_w1 <- sweep(log_w1, c(1, 3), apply(log_w1, c(1, 3), max))
    wp0 <- exp(log_w0) * post_prob
    wp1 <- exp(log_w1) * post_prob
    
    hzrd0_base <- exp(sweep(Y0, c(1,3), theta0, FUN = '+'))
    hzrd1_base <- exp(sweep(Y1, c(1,3), theta1, FUN = '+'))
    hzrd0 <- sweep(hzrd0_base, c(1, 3), t^(exp(theta0) - 1), FUN = '*')
    hzrd1 <- sweep(hzrd1_base, c(1, 3), t^(exp(theta1) - 1), FUN = '*')
    sum_hzrd0 <- apply(hzrd0 * wp0, c(1, 3), sum)
    sum_hzrd1 <- apply(hzrd1 * wp1, c(1, 3), sum)
    sum_wp0 <- apply(wp0, c(1, 3), sum)
    sum_wp1 <- apply(wp1, c(1, 3), sum)
    mean_hzrd0 <- sum_hzrd0 / sum_wp0
    mean_hzrd1 <- sum_hzrd1 / sum_wp1
    hazard_ratio <- mean_hzrd0 / mean_hzrd1
    return (list(
      Hazard_0 = mean_hzrd0,
      Hazard_1 = mean_hzrd1,
      Hazard_Ratio = hazard_ratio
    ))
  }
  return (structure(list(
    Prob = post_prob_raw,
    func = f
    ), class = "post_survival"))
}

plot.post_survival <- function(result, t_range = c(0, 10), t_count = 100) {
  get_long_df <- function(var, data) {
    df <- tidyr::pivot_longer(as.data.frame(t(data[[var]])), everything())
    df$var <- var
    return (df)
  }
  result$Prob <- apply(result$Prob, c(1, 3), mean)
  long_prob <- get_long_df("Prob", result)
  t_ticks <- seq(max(t_range[1], t_range[2] / 100), t_range[2], length.out = t_count)
  get_df <- function(t) {
    r <- result$func(t)
    h0 <- get_long_df("Hazard_0", r)
    h1 <- get_long_df("Hazard_1", r)
    hr <- get_long_df("Hazard_Ratio", r)
    return (dplyr::bind_rows(h0, h1, hr))
  }
  all_df <- dplyr::bind_rows(pbapply::pblapply(t_ticks, function(t) bind_cols(get_df(t), t = t)))
  
  p0 <- ggplot(long_prob) + 
    geom_density(aes(x = value, color = name)) +
    geom_histogram(aes(x = value, y = after_stat(density), fill = name), alpha = 0.3) +
    guides(color = "none", fill = "none") +
    ggtitle("Probability") + xlab("probability")
  df_1 <- all_df %>% filter(var != "Hazard_Ratio") %>%
    group_by(name, var, t) %>%
    summarize(mean = mean(value), lwr = quantile(value, 0.025), upr = quantile(value, 0.975))
  p1 <- ggplot(df_1) + 
    geom_line(aes(x = t, y = mean, color = var)) +
    geom_ribbon(aes(x = t, ymin = lwr, ymax = upr, fill = var), alpha = 0.3) + 
    facet_wrap(~name) + ggtitle("Hazard") + xlab("time") + ylab("hazard")
  
  df_2 <- all_df %>% filter(var == "Hazard_Ratio") %>%
    group_by(name, var, t) %>%
    summarize(mean = mean(value), lwr = quantile(value, 0.025), upr = quantile(value, 0.975))
  p2 <- ggplot(df_2) + 
    geom_line(aes(x = t, y = mean, color = name)) +
    geom_ribbon(aes(x = t, ymin = lwr, ymax = upr, fill = name), alpha = 0.3) +
    ggtitle("Hazard Ratio") + xlab("time") + ylab("hazard ratio")
  (p0 + p2) / p1
}

plot.post_outcome <- function(result) {
  get_long_df <- function(var) {
    df <- tidyr::pivot_longer(as.data.frame(t(result[[var]])), everything())
    df$var <- var
    return (df)
  }
  result$Prob <- apply(result$Prob, c(1, 3), mean)
  long_prob <- get_long_df("Prob")
  long_O1 <- get_long_df("Overall_1")
  long_O0 <- get_long_df("Overall_0")
  long_CE <- get_long_df("Causal_Effect")
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

print.PStrata <- function(res){
  if (res$PSobject$Y.family$family == "survival")
    cat("Survival Data\n")
  cat("Posterior estimate of the parameters:\n")
  mat <- rstan::summary(res$stanfit)$summary
  #param_names <- sapply(obj$PSobject$parameter_list, function(x) x$name)
  #rownames(mat) <- c(param_names, "__lp")
  print(mat)
  cat('\n')
  cat("Estimated Proportion from Each Stratum:\n")
  prob_samples <- apply(res$outcome$Prob, c(1, 3), mean)
  prob_summary <- apply(prob_samples, 1, function(x) {
    c(mean = mean(x), sd = sd(x), lwr = quantile(x, 0.025), upr = quantile(x, 0.975))
  })
  print(prob_summary)
  cat('\n')
  if (res$PSobject$Y.family$family != "survival") {
    cat("Causal Effect of Principal Strata:\n")
    outcome_summary <- apply(res$outcome$Causal_Effect, 1, function(x) {
      c(mean = mean(x), sd = sd(x), lwr = quantile(x, 0.025), upr = quantile(x, 0.975))
    })
    print(outcome_summary)
    cat('\n')
  }
  else {
    cat("Hazard Ratio of Principal Strata at t = 1:\n")
    outcome_summary <- apply(res$outcome$func(1)$Hazard_Ratio, 1, function(x) {
      c(mean = mean(x, na.rm = T), sd = sd(x, na.rm = T), 
        lwr = quantile(x, 0.025, na.rm = T), upr = quantile(x, 0.975, na.rm = T))
    })
    print(outcome_summary)
    cat('\n')
  }
}

plot.PStrata <- function(res, ...){
  plot(res$outcome, ...)
}
