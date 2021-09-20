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
  df <- list(
    N = nrow(data), 
    Z = dplyr::pull(data, obj$variables$treatment),
    D = dplyr::pull(data, obj$variables$intervention),
    Y = dplyr::pull(data, obj$variables$outcome),
    X = dplyr::select(data, as.character(obj$symbol_list))
  )
  if (!is.null(obj$variables$censor))
    df$C <- dplyr::pull(data, obj$variables$censor)
  stanfit <- rstan::stan(paste0(model_name, ".stan"), data = df, 
                         ...
  )
  pso_code <- paste(readLines(paste0(model_name, '.pso')), collapse = '\n')
  stan_code <- paste(readLines(paste0(model_name, '.stan')), collapse = '\n')
  posterior_samples <- do.call(rbind, rstan::extract(stanfit))
  
  avg <- function(par) {
    mat <- as.matrix(df$X)
    single_outcome <- function(stratum) {
      mean1 <- apply(mat, 1, function(x) obj$outcome_model_list[[stratum]][[2]]$eval(x, par))
      mean0 <- apply(mat, 1, function(x) obj$outcome_model_list[[stratum]][[1]]$eval(x, par))
      prob <- apply(mat, 1, function(x) obj$stratum_model_list[[stratum]]$eval(x, par))
      prob <- prob - max(prob)
      if (Y.family == "survival"){
        mean1 <- exp(mean1)
        mean0 <- exp(mean0)
      }
      res1 <- sum(mean1 * exp(prob)) / sum(exp(prob))
      res0 <- sum(mean0 * exp(prob)) / sum(exp(prob))
      return (c(res1, res0))
    }
    resmat <- matrix(nrow = 4, ncol = length(obj$strata))
    dimnames(resmat) <- list(c('Z = 0', 'Z = 1', 'ATE', 'Prob'), obj$strata)
    for (stratum in obj$strata){
      resmat[1:2, stratum] <- single_outcome(stratum)
    }
    
    if (Y.family == "survival")
      resmat[3, ] <- resmat[1, ] / resmat[2, ] # ? should do this?
    else
      resmat[3, ] <- resmat[1, ] - resmat[2, ]
    
    # probability
    get_prob <- function(x) {
      log_probs <- sapply(obj$strata, function(s) obj$stratum_model_list[[s]]$eval(x, par))
      log_probs <- log_probs - max(log_probs)
      return (exp(log_probs) / sum(exp(log_probs)))
    }
    resmat[4, ] <- rowMeans(apply(mat, 1, get_prob))
    
    return (resmat)
  }
  
  causal_effects <- pbapply::pbapply(posterior_samples, 2, avg)
  mean_effect <- apply(causal_effects, 1, mean)
  se_effect <- apply(causal_effects, 1, sd)
  effect_table <- matrix(nrow = 8, ncol = length(obj$strata))
  dimnames(effect_table) <- list(
    c("Prob (mean)", "Prob (sd)", "Z=1 (mean)", "Z=1 (sd)", "Z=0 (mean)", "Z=0 (sd)", "ATE (mean)", "ATE (sd)"), obj$strata
  )
  effect_table[c(3, 5, 7, 1), ] <- matrix(mean_effect, nrow = 4, byrow = F)
  effect_table[c(4, 6, 8, 2), ] <- matrix(se_effect, nrow = 4, byrow = F)
  
  param_names <- c(sapply(obj$parameter_list, function(x) x$name), "__lp")
  rownames(posterior_samples) <- param_names
  
  res <- list(
    data = df,
    stanfit = stanfit,
    PSobject = obj,
    pso_code = pso_code,
    stan_code = stan_code,
    post_samples = posterior_samples,
    causal_effects = causal_effects,
    effect_table = effect_table
  )
  class(res) <- "PStrata"
  return (res)
}

print.PStrata <- function(obj){
  cat("Posterior estimate of the parameters:\n")
  mat <- rstan::summary(obj$stanfit)$summary
  param_names <- sapply(obj$PSobject$parameter_list, function(x) x$name)
  rownames(mat) <- c(param_names, "__lp")
  print(mat)
  cat('\n')
  cat("Estimated Proportion from Each Stratum:\n")
  n_strata <- length(obj$PSobject$strata)
  prop_S <- obj$causal_effects[4 * (1:n_strata), ]
  mat1 <- t(apply(prop_S, 1, function(x) c(mean(x), sd(x),
                                           quantile(x, c(0.025, 0.975)))))
  colnames(mat1) <- c("mean", "sd", "lwr", "upr")
  rownames(mat1) <- obj$PSobject$strata
  print(mat1)
  
  cat('\n')
  cat("Causal Effect of Principal Strata:\n")
  causal_Y <- obj$causal_effects[4 * (1:n_strata) - 1, ]
  mat2 <- t(apply(causal_Y, 1, function(x) c(mean(x), sd(x),
                                             quantile(x, c(0.025, 0.975)))))
  colnames(mat2) <- c("mean", "sd", "lwr", "upr")
  rownames(mat2) <- obj$PSobject$strata
  print(mat2)
}

plot.PStrata <- function(obj){
  data_raw <- data.frame(t(obj$post_samples))
  data_long <- do.call(dplyr::bind_rows, lapply(names(data_raw), 
                                                function(col) data.frame(
                                                  x = density(data_raw[, col])$x, 
                                                  y = density(data_raw[, col])$y, 
                                                  fill = density(data_raw[, col])$x > quantile(data_raw[, col], .025) &
                                                    density(data_raw[, col])$x < quantile(data_raw[, col], .975),
                                                  name = col)))
  data_long$name <- factor(data_long$name, levels = colnames(data_raw))
  ggplot2::ggplot(data_long) + 
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymax = y, ymin = 0), alpha = 0.9, fill = "springgreen2") + 
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymax = y, ymin = 0), alpha = 0.9, fill = "springgreen4", data = data_long[data_long$fill, ]) +
    ggplot2::facet_wrap(~name, scales = "free")
}
