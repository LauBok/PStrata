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
    df$C = dplyr::pull(data, obj$variables$censor)
  stanfit <- rstan::stan(paste0(model_name, ".stan"), data = df, 
                         ...
  )
  pso_code <- paste(readLines(paste0(model_name, '.pso')), collapse = '\n')
  stan_code <- paste(readLines(paste0(model_name, '.stan')), collapse = '\n')
  return (stanfit)
  
  posterior_samples <- do.call(rbind, rstan::extract(stanfit))
  get_mean_difference_individual <- function(i, j){
    cur_effect <- matrix(0, nrow = 3, ncol = length(obj$strata), 
                         dimnames = list(c("effect", "weight", "bool"), obj$strata))
    for (stratum in obj$strata){
      if (as.numeric(substring(stratum, df$Z[i] + 1, df$Z[i] + 1)) == df$D[i]){
        # this stratum counts
        cur_effect["bool", stratum] = 1
        # calculate the probability
        log_prob <- obj$S_list[[stratum]]$fun(df$X[i,], posterior_samples[, j])
        cur_effect["weight", stratum] <- log_prob
        # calculate the causal effect
        tmp_obj <- obj$Y_list[[stratum]]
        out_1 <- tmp_obj[[2]]$fun(df$X[i,], posterior_samples[, j])
        out_0 <- tmp_obj[[1]]$fun(df$X[i,], posterior_samples[, j])
        out_1_inv_link <- tmp_obj[[2]]$invlink(out_1)
        out_0_inv_link <- tmp_obj[[1]]$invlink(out_0)
        cur_effect["effect", stratum] <-out_1_inv_link - out_0_inv_link
      }
    }
    cur_effect["weight", ] = cur_effect["weight", ] - max(cur_effect["weight", ])
    cur_effect["weight", ] = exp(cur_effect["weight", ]) * cur_effect["bool", ]
    cur_effect["weight", ] = cur_effect["weight", ] / sum(cur_effect["weight", ])
    cur_effect["effect", ] = cur_effect["effect", ] * cur_effect["weight", ]
    return (cur_effect[c("effect", "weight"), ])
  }
  
  get_mean_difference <- function(j){
    causal_effect <- matrix(0, nrow = 2, ncol = length(obj$strata), 
                            dimnames = list(c("effect", "weight"), obj$strata))
    for (i in 1:df$N){
      causal_effect <- causal_effect + get_mean_difference_individual(i, j)
    }
    return (c(causal_effect["weight", ], causal_effect["effect", ] / causal_effect["weight", ]))
  }
  
  causal_effects <- pbapply::pbsapply(1:ncol(posterior_samples), get_mean_difference)
  
  param_names <- c(sapply(obj$param_list, function(x) x$name), "__lp")
  rownames(posterior_samples) <- param_names
  
  res <- list(
    data = df,
    stanfit = stanfit,
    PSobject = obj,
    pso_code = pso_code,
    stan_code = stan_code,
    post_samples = posterior_samples,
    causal_effect = causal_effects
  )
  class(res) <- "PStrata"
  #return (res)
}

print_future.PStrata <- function(obj){
  cat("Posterior estimate of the parameters:\n")
  mat <- rstan::summary(obj$stanfit)$summary
  param_names <- sapply(obj$PSobject$param_list, function(x) x$name)
  rownames(mat) <- c(param_names, "__lp")
  print(mat)
  cat('\n')
  cat("Estimated Proportion from Each Stratum:\n")
  prop_S <- obj$causal_effect[1 : (nrow(obj$causal_effect) %/% 2), ]
  prop_S <- apply(prop_S, 2, function(x) x / sum(x))
  mat1 <- t(apply(prop_S, 1, function(x) c(mean(x), sd(x),
                                           quantile(x, c(0.025, 0.975)))))
  colnames(mat1) <- c("mean", "sd", "lwr", "upr")
  print(mat1)
  
  cat('\n')
  cat("Causal Effect of Principal Strata:\n")
  causal_effect_Y <- obj$causal_effect[(1 + nrow(obj$causal_effect) %/% 2) : nrow(obj$causal_effect), ]
  mat2 <- t(apply(causal_effect_Y, 1, function(x) c(mean(x), sd(x),
                                                    quantile(x, c(0.025, 0.975)))))
  colnames(mat2) <- c("mean", "sd", "lwr", "upr")
  print(mat2)
}

plot_future.PStrata <- function(obj){
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
