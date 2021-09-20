### gaussian
build_gaussian <- function(link, param_start_num){
  .param_list <- list(list(
    type = "positive", prior_type = "prior_sigma", name = "Sigma"
  ))
  .param <- str2lang(paste0('`@', param_start_num, '`'))
  .mean_func <- function(core){
    if (link == "identity")
      return (core)
    if (link == "log")
      return (as.call(list(quote(exp), core)))
    if (link == "inverse")
      return (as.call(list(quote(`/`), 1, core)))
  }
  .model_func <- function(core){
    return (as.call(list(
      quote(normal_lpdf), .mean_func(core), .param
    )))
  }
  return (list(
    param_list = .param_list,
    mean_func = .mean_func,
    model_func = .model_func
  ))
}

### binomial
build_binomial <- function(link, param_start_num){
  .param_list <- list()
  .mean_func <- function(core){
    if (link == "logit")
      return (as.call(list(quote(inv_logit), core)))
    if (link == "probit")
      return (as.call(list(quote(inv_Phi), core)))
    if (link == "cauchit"){
      .atan <- as.call(list(quote(atan), core))
      .pi <- as.call(list(quote(pi)))
      .div <- as.call(list(quote(`/`), .atan, .pi))
      return (as.call(list(quote(`+`), .5, .div)))
    }
    if (link == "log")
      return (as.call(list(quote(exp), core)))
    if (link == "cloglog")
      return (as.call(list(quote(inv_cloglog), core)))
  }
  .model_func <- function(core){
    return (as.call(list(
      quote(bernoulli_lpmf), .mean_func(core)
    )))
  }
  return (list(
    param_list = .param_list,
    mean_func = .mean_func,
    model_func = .model_func
  ))
}

### Gamma
build_Gamma <- function(link, param_start_num){
  .param_list <- list(list(
    type = "positive", prior_type = "prior_alpha", name = "Alpha"
  ))
  .param <- str2lang(paste0('`@', param_start_num, '`'))
  .mean_func <- function(core){
    if (link == "identity")
      return (as.call(list(quote(`/`), .param, core)))
    if (link == "log")
      return (as.call(list(quote(`/`), .param, 
                           as.call(list(quote(exp), core)))))
    if (link == "inverse")
      return (as.call(list(quote(`*`), .param, core)))
  }
  .model_func <- function(core){
    return (as.call(list(
      quote(gamma_lpdf), .param, .mean_func(core)
    )))
  }
  return (list(
    param_list = .param_list,
    mean_func = .mean_func,
    model_func = .model_func
  ))
}

### poisson
build_poisson <- function(link, param_start_num){
  .param_list <- list()
  .mean_func <- function(core){
    if (link == "identity")
      return (core)
    if (link == "log")
      return (as.call(quote(exp), core))
    if (link == "sqrt")
      return (as.call(quote(`^`), core, 2))
  }
  .model_func <- function(core){
    return (as.call(list(
      quote(poisson_lpmf), .mean_func(core)
    )))
  }
  return (list(
    param_list = .param_list,
    mean_func = .mean_func,
    model_func = .model_func
  ))
}

### inverse.gaussian
build_inv_guassian <- function(link, param_start_num){
  .param_list <- list(list(
    type = "positive", prior_type = "prior_lambda", name = "Lambda"
  ))
  .param <- str2lang(paste0('`@', param_start_num, '`'))
  .mean_func <- function(core){
    if (link == "identity")
      return (core)
    if (link == "log")
      return (as.call(list(quote(exp), core)))
    if (link == "inverse")
      return (as.call(list(quote(`/`), 1, core)))
    if (link == "1/mu^2")
      return (as.call(list(quote(`^`), core, -2)))
  }
  .model_func <- function(core){
    return (as.call(list(
      quote(inv_gaussian_lpdf), .mean_func(core), .param
    )))
  }
  return (list(
    param_list = .param_list,
    mean_func = .mean_func,
    model_func = .model_func
  ))
}

survival <- function() {
  return (list(
    family = "survival"
  ))
}

### survival
build_survival <- function(link, param_start_num){
  .param_list <- list(
    list(type = "real", prior_type = "prior_theta", name = "Theta1"),
    list(type = "real", prior_type = "prior_theta", name = "Theta2")
  )
  .param1 <- str2lang(paste0('`@', param_start_num, '`'))
  .param2 <- str2lang(paste0('`@', param_start_num + 1, '`'))
  .mean_func <- function(core){
    return (
      as.call(list(quote(`+`), as.call(list(quote(`+`), .param1, .param2)), core))
    )
  }
  .model_func <- function(core){
    return (as.call(list(
      quote(survival_lpdf), core, .param1, .param2, quote(.)
    )))
  }
  return (list(
    param_list = .param_list,
    mean_func = .mean_func,
    model_func = .model_func
  ))
}

build_model <- function(family, param_start_num = 1L){
  if (family$family == "gaussian")
    res <- build_gaussian(family$link, param_start_num)
  else if (family$family == "binomial")
    res <- build_binomial(family$link, param_start_num)
  else if (family$family == "Gamma")
    res <- build_Gamma(family$link, param_start_num)
  else if (family$family == "poisson")
    res <- build_poisson(family$link, param_start_num)
  else if (family$family == "inverse.gaussian")
    res <- build_inv_gaussian(family$link, param_start_num)
  else if (family$family == "survival")
    res <- build_survival(family$link, param_start_num)
  
  return (list(
    params = res$param_list,
    param_cnt = length(res$param_list),
    mean = res$mean_func,
    model = res$model_func
  ))
}

model_template <- function(model, family, ER, stratum){
  # OUTPUT: generated parameter list
  # OUTPUT: generated model list
  .func <- function(start_num = 1L){
    param_list <- list()
    model_list <- list()
    
    # model params
    # this is the set of parameters related to the model
    # we separate this part out because it might need to be reused
    # depending whether or not ER is assumed
    model_result <- manipulate_formula(model)
    
    for (i in 1:ifelse(ER, 1, 2)){
      # generate a set of parameters for the model
      # only generate a second set if ER is not assumed
      core <- model_result$evaluate(start_num)
      new_params <- model_result$parameters
      if (length(new_params) != 0){
        for (j in 1:(length(new_params))){
          new_params[[j]]$name <- paste(
            'Y', stratum, i - 1, new_params[[j]]$name, sep = "_"
          )
        }
      }
      param_list <- c(param_list, new_params)
      start_num <- start_num + model_result$num_of_parameters
      model.tmp <- build_model(family, start_num)
      new_params <- model.tmp$params
      if (length(new_params) != 0){
        for (j in 1:(length(new_params))){
          new_params[[j]]$name <- paste(
            'Y', stratum, i - 1, new_params[[j]]$name, sep = "_"
          )
        }
      }
      param_list <- c(param_list, new_params)
      model_list[[i]] <- list(
        mean = model.tmp$mean(core),
        model = model.tmp$model(core)
      )
      start_num <- start_num + model.tmp$param_cnt
    }
    
    if (ER)
      model_list[[2]] <- model_list[[1]]
    
    return (list(
      param_list = param_list,
      model_list = model_list,
      start = start_num
    ))
  }
  
  return (.func)
}