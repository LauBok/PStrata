get_data <- function(PSobject, data) {
  df <- list(
    N = nrow(data), 
    Z = dplyr::pull(data, PSobject$variables$treatment),
    D = dplyr::pull(data, PSobject$variables$intervention),
    Y = dplyr::pull(data, PSobject$variables$outcome)
  )
  if (!is.na(PSobject$variables$censor))
    df$C <- dplyr::pull(data, PSobject$variables$censor)
  if (PSobject$S.dim){
    df$XS <- model.matrix(update.formula(PSobject$S.model, . ~ . - 1), data)
  }
  if (PSobject$Y.dim){
    df$XY <- model.matrix(update.formula(PSobject$Y.model, . ~ . - 1), data)
  }
  return (df)
}

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

get_posterior_samples <- function(stanfit, obj) {
  posterior_samples <- do.call(cbind, rstan::extract(stanfit))
  param_names <- c()
  for (param in obj$parameter_list) {
    param_names <- c(param_names, generate_names(param$class, param$group, param$name, param$dim))
  }
  param_names <- c(param_names, '__lp')
  colnames(posterior_samples) <- param_names
  return (posterior_samples)
}

PStrata <- function(S.formula, Y.formula, Y.family, data, monotonicity = "default", ER = c(), trunc = FALSE, 
               prior_intercept = prior_uniform(),
               prior_coefficient = prior_normal(),
               prior_sigma = prior_inv_gamma(),
               prior_alpha = prior_inv_gamma(),
               prior_lambda = prior_inv_gamma(),
               prior_theta = prior_normal(),
               model_name = "unnamed", ...){
  obj <- PSObject(
    S.formula, Y.formula, Y.family, monotonicity, ER,
    prior_intercept, prior_coefficient, prior_sigma,
    prior_alpha, prior_lambda, prior_theta
  )
  write.pso(obj, paste0(model_name, ".pso"))
  to_stan(model_name)
  dataset <- df <- get_data(obj, data)
  stanfit <- rstan::stan(paste0(model_name, ".stan"), data = df, 
                         ...
  )
  pso_code <- paste(readLines(paste0(model_name, '.pso')), collapse = '\n')
  stan_code <- paste(readLines(paste0(model_name, '.stan')), collapse = '\n')
  post_samples <- get_posterior_samples(stanfit, obj)
  post_probability_raw <- post_prob_raw(obj, post_samples, df)
  post_probability <- post_prob(obj, post_samples, df)
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

post_prob_raw <- function(obj, post_samples, df){
  results <- array(
    0, 
    dim = c(length(obj$strata), length(df$Y), nrow(post_samples)),
    dimnames = list(obj$strata, NULL, NULL)
  ) # log_probability
  names <- colnames(post_samples)
  for (stratum in obj$strata[-1]) {
    value <- 0
    if (generate_names('S', stratum, 'Intrcpt') %in% names)
      intercept <- post_samples[, generate_names('S', stratum, 'Intrcpt'), drop = F]
    else
      intercept <- 0
    value <- value + rep(1, each = length(df$Y)) %*% t(intercept)
    if (obj$S.dim > 0){
      pars <- post_samples[, generate_names('S', stratum, 'Coef', obj$S.dim)]
      value <- value + df$XS %*% t(pars)
    }
    results[stratum, , ] <- value
    
  }
  max_log_prob <- apply(results, c(2, 3), max)
  results <- exp(sweep(results, c(2, 3), max_log_prob))
  sum_prob <- apply(results, c(2, 3), sum)
  results <- sweep(results, c(2, 3), sum_prob, FUN = '/')
  return (results)
}

post_prob <- function(obj, post_samples, df) {
  results <- post_prob_raw(obj, post_samples, df)
  consistency <- array(1, dim = dim(results),
                       dimnames = list(obj$strata, NULL, NULL))
  for (stratum in obj$strata[-1]) {
    # consistency
    consistency[stratum, , ] <- ifelse(substring(stratum, df$Z+1, df$Z+1) == df$D, 1, 0)
  }
  results <- results * consistency
  sum_prob <- apply(results, c(2, 3), sum)
  results <- sweep(results, c(2, 3), sum_prob, FUN = '/')
  return (results)
}

post_Y <- function(obj, post_samples, df) {
  results <- array(
    0, 
    dim = c(length(obj$strata), 2, length(df$Y), nrow(post_samples)),
    dimnames = list(obj$strata, c("0", "1"), NULL, NULL)
  ) # log_probability
  theta <- array(0, dim = c(length(obj$strata), 2, nrow(post_samples)),
                 dimnames = list(obj$strata, c("0", "1"), NULL))
  names <- colnames(post_samples)
  for (stratum in obj$strata) {
    for (z in c("0", "1")){
      if (stratum %in% obj$ER)
        zz = "0"
      else 
        zz = z
      value <- 0
      if (generate_names('Y', c(stratum, zz), 'Intrcpt') %in% names)
        intercept <- post_samples[, generate_names('Y', c(stratum, zz), 'Intrcpt'), drop = F]
      else
        intercept <- 0
      value <- value + rep(1, each = length(df$Y)) %*% t(intercept)
      if (obj$S.dim > 0){
        pars <- post_samples[, generate_names('Y', c(stratum, zz), 'Coef', obj$Y.dim)]
        value <- value + df$XY %*% t(pars)
      }
      if (obj$Y.family$family == "survival") {
        theta[stratum, z, ] <- post_samples[, generate_names('Y', c(stratum, zz), 'Theta')]
      }
      
      results[stratum, z, , ] <- value
    }
  }
  if (obj$Y.family$family != "survival"){
    return (obj$Y.family$linkinv(results))
  }
  else {
    return (list(mean = results, theta = theta))
  }
}

post_stratum_samples <- function(obj, post_samples, df) {
  prob <- post_prob(obj, post_samples, df)
  apply(prob, c(2,3), function(p) (1:dim(prob)[1]) %*% rmultinom(1,1,p))
}

plot_prob_prob_one <- function(num_draw, obj, post_samples, df) {
  post_stratum <- post_stratum_samples(obj, post_samples, df)[, num_draw]
  prob <- post_prob(obj, post_samples, df)[,, num_draw]
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
