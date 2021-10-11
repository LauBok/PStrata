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
    DimS = obj$S.dim,
    DimY = obj$Y.dim,
    XS = model.matrix(update.formula(obj$S.model, . ~ . - 1), data),
    XY = model.matrix(update.formula(obj$Y.model, . ~ . - 1), data)
  )
  if (!is.null(obj$variables$censor))
    df$C <- dplyr::pull(data, obj$variables$censor)
  stanfit <- rstan::stan(paste0(model_name, ".stan"), data = df, 
                         ...
  )
  pso_code <- paste(readLines(paste0(model_name, '.pso')), collapse = '\n')
  stan_code <- paste(readLines(paste0(model_name, '.stan')), collapse = '\n')
  posterior_samples <- do.call(rbind, rstan::extract(stanfit))
  
  res <- list(
    data = df,
    stanfit = stanfit,
    PSobject = obj,
    pso_code = pso_code,
    stan_code = stan_code,
    post_samples = posterior_samples,
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
