PS <- function(S.formula, Y.formula, Y.family, data, monotonicity = "default", ER = c(), trunc = FALSE, ...){
  obj <- PSobject(S.formula, Y.formula, Y.family, monotonicity, ER, trunc)
  write.pso(obj, "model.pso")
  main()
  df <- list(
    N = nrow(data), 
    Z = data[, obj$treatment_var],
    D = data[, obj$intervention_var],
    Y = data[, obj$outcome_var],
    X = as.matrix(data[, obj$symbol_list, drop = F])
  )
  df_X <- df$X
  stanfit <- rstan::stan("model.stan", data = df, 
      ...
  )
  pso_code <- paste(readLines('model.pso'), collapse = '\n')
  stan_code <- paste(readLines('model.stan'), collapse = '\n')
  #file.remove('model.pso')
  #file.remove('model.stan')
  
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
    return (causal_effect["effect", ] / causal_effect["weight", ])
  }
  
  causal_effects <- pbapply::pbsapply(1:ncol(posterior_samples), get_mean_difference)
  
  
  res <- list(
    stanfit = stanfit,
    PSobject = obj,
    pso_code = pso_code,
    stan_code = stan_code,
    post_samples = posterior_samples,
    causal_effect = causal_effects
  )
  class(res) <- "PSstan"
  return (res)
}

print.PSstan <- function(obj){
  cat("Posterior estimate of the parameters:\n")
  mat <- rstan::summary(obj$stanfit)$summary
  param_names <- sapply(obj$PSobject$param_list, function(x) x$name)
  rownames(mat) <- c(param_names, "__lp")
  print(mat)
  cat('\n')
  cat("Causal Effect of Principal Strata:\n")
  mat2 <- t(apply(obj$causal_effect, 1, function(x) c(mean(x), sd(x),
                                                      quantile(x, c(0.025, 0.975)))))
  colnames(mat2) <- c("mean", "sd", "lwr", "upr")
  print(mat2)
}
