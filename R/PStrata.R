#' Principal strata analysis for data with post-randomization intervention
#' 
#' description
#' 
#' @export PStrata
#' @useDynLib PStrata
#' @param S.formula a formula for the principal strata model, in the form of \code{Z + D1 + D2 ~ X1 + X2}, where \code{Z} is the treatment assignment and \code{D1} and \code{D2} are the post-randomization variables.
#' @param Y.formula a formula for the outcome model, in the form of \code{Y ~ X1 + X2}, where \code{Y} is the outcome variable
#' @param Y.family a family object, indicating the family of the model as used in \code{lm()}.
#' @param data a data frame object, the data used to do inference on.
#' @param strata a vector of integers, each integer referring to a stratum that is included in the model
#' @param ER a vector of integers, indicating the strata on which the exclusion restriction assumption is assumed.
#' @param prior_... the prior distribution for the corresponding parameters.
#' @return A text file, readable by \code{to_stan()}.

PStrata <- function(
    S.formula, Y.formula, Y.family, data, strata, ER,
    prior_intercept = prior_flat(),
    prior_coefficient = prior_normal(),
    prior_sigma = prior_inv_gamma(),
    prior_alpha = prior_inv_gamma(),
    prior_lambda = prior_inv_gamma(),
    prior_theta = prior_normal(),
    survival.time.points = 50,
    filename = "unnamed",
    ...
) {
  PSobject <- PSObject(S.formula, Y.formula, Y.family, data, strata, ER,
                       prior_intercept,
                       prior_coefficient,
                       prior_sigma,
                       prior_alpha,
                       prior_lambda,
                       prior_theta,
                       survival.time.points,
                       filename)
  write.auxillary.files()
  to_stan(filename)
  post_samples <- PSSample(filename = paste0(filename, ".stan"),
           data = PSobject$stan_data,
           ...)
  res <- list(
    PSobject = PSobject,
    post_samples = post_samples
  )
  if (Y.family$family == "survival_Cox" || Y.family$family == "survival_AFT")
    class(res) <- "PStrata_survival"
  else
    class(res) <- "PStrata"
  return (res)
}

#' @method  summary PStrata
#' @export summary.PStrata
summary.PStrata <- function(pstrata) {
  samples <- pstrata$post_samples
  strata_prob_samples <- rstan::extract(samples)$`strata_prob`
  mean_effect_samples <- rstan::extract(samples)$`mean_effect`
  treatment_effect_samples <- 0 * strata_prob_samples
  strata_prob_summary <- t(apply(strata_prob_samples, 2, function(x) c(
    mean = mean(x), sd = sd(x),
    `2.5%` = unname(quantile(x, 0.025)), `25%` = unname(quantile(x, 0.25)),
    `median` = unname(quantile(x, 0.5)),
    `75%` = unname(quantile(x, 0.75)), `97.5%` = unname(quantile(x, 0.975))
  )))
  mean_effect_summary <- t(apply(mean_effect_samples, 2, function(x) c(
    mean = mean(x), sd = sd(x),
    `2.5%` = unname(quantile(x, 0.025)), `25%` = unname(quantile(x, 0.25)),
    `median` = unname(quantile(x, 0.5)),
    `75%` = unname(quantile(x, 0.75)), `97.5%` = unname(quantile(x, 0.975))
  )))
  meta_data <- pstrata$PSobject$meta_data
  rownames(strata_prob_summary) <- paste("S = ", meta_data$strata, sep = '')
  group_names <- rep(NA, nrow(mean_effect_summary))
  for (i in 1:nrow(meta_data$SZDG)){
    cur_g = meta_data$SZDG[i, 4] + 1
    cur_s = meta_data$SZDG[i, 1] + 1
    cur_z = meta_data$SZDG[i, 2]
    if (!is.na(group_names[cur_g])) {
      group_names[cur_g] <- 
        paste0("S = ", meta_data$strata[cur_s], ", Z = *")
    }
    else {
      group_names[cur_g] <- 
        paste0("S = ", meta_data$strata[cur_s], ", Z = ", cur_z)
    }
  }
  rownames(mean_effect_summary) <- group_names
  
  for (j in 1:ncol(treatment_effect_samples)) {
    g1 <- meta_data$SZDG[2 * j, 4] + 1
    g0 <- meta_data$SZDG[2 * j - 1, 4] + 1
    treatment_effect_samples[, j] <- 
      mean_effect_samples[, g1] - mean_effect_samples[, g0]
  }
  treatment_effect_summary <- t(apply(treatment_effect_samples, 2, function(x) c(
    mean = mean(x), sd = sd(x),
    `2.5%` = unname(quantile(x, 0.025)), `25%` = unname(quantile(x, 0.25)),
    `median` = unname(quantile(x, 0.5)),
    `75%` = unname(quantile(x, 0.75)), `97.5%` = unname(quantile(x, 0.975))
  )))
  rownames(treatment_effect_summary) <- paste("S = ", meta_data$strata, sep = '')
  
  return (list(strata_prob = strata_prob_summary,
               mean_effect = mean_effect_summary,
               treatment_effect = treatment_effect_summary))
}

#' @method  plot PStrata
#' @export plot.PStrata
plot.PStrata <- function(pstrata, 
                         type = c("strata_prob", "mean_effect", "treatment_effect")) {
  samples <- pstrata$post_samples
  meta_data <- pstrata$PSobject$meta_data
  strata_prob_samples <- rstan::extract(samples)$`strata_prob`
  colnames(strata_prob_samples) <- meta_data$strata
  strata_prob_plot <- tidyr::pivot_longer(as.data.frame(strata_prob_samples), 
                                          everything(), names_to = 'Stratum',
                                          values_to = 'Probability')
  mean_effect_samples <- rstan::extract(samples)$`mean_effect`
  group_names <- rep(NA, ncol(mean_effect_samples))
  for (i in 1:nrow(meta_data$SZDG)){
    cur_g = meta_data$SZDG[i, 4] + 1
    cur_s = meta_data$SZDG[i, 1] + 1
    cur_z = meta_data$SZDG[i, 2]
    if (!is.na(group_names[cur_g])) {
      group_names[cur_g] <- 
        paste0("S = ", meta_data$strata[cur_s], ", Z = *")
    }
    else {
      group_names[cur_g] <- 
        paste0("S = ", meta_data$strata[cur_s], ", Z = ", cur_z)
    }
  }
  colnames(mean_effect_samples) <- group_names
  mean_effect_plot <- tidyr::pivot_longer(as.data.frame(mean_effect_samples), 
                                          everything(), names_to = 'Group',
                                          values_to = 'Mean Effect')
  treatment_effect_samples <- 0 * strata_prob_samples
  for (j in 1:ncol(treatment_effect_samples)) {
    g1 <- meta_data$SZDG[2 * j, 4] + 1
    g0 <- meta_data$SZDG[2 * j - 1, 4] + 1
    treatment_effect_samples[, j] <- 
      mean_effect_samples[, g1] - mean_effect_samples[, g0]
  }
  colnames(treatment_effect_samples) <- meta_data$strata
  treatment_effect_plot <- tidyr::pivot_longer(as.data.frame(treatment_effect_samples), 
                                          everything(), names_to = 'Stratum',
                                          values_to = 'Treatment Effect')
  
  Gplot <- NULL
  if (type == "strata_prob" || type == 1) {
    Gplot <- ggplot2::ggplot(strata_prob_plot, ggplot2::aes(Probability)) +
      ggplot2::geom_density(ggplot2::aes(fill = Stratum), alpha = 0.4)
  }
  else if (type == "mean_effect" || type == 2) {
    Gplot <- ggplot2::ggplot(mean_effect_plot, ggplot2::aes(`Mean Effect`)) +
      ggplot2::geom_density(ggplot2::aes(fill = Group), alpha = 0.4)
  }
  else if (type == "treatment_effect" || type == 3) {
    Gplot <- ggplot2::ggplot(treatment_effect_plot, ggplot2::aes(`Treatment Effect`)) +
      ggplot2::geom_histogram(ggplot2::aes(fill = Stratum, y = ..density..), alpha = 0.4)
  }
  
  return (Gplot)
}

#' @export contrast
contrast <- function(pstrata) {
  samples <- pstrata$post_samples
  strata_prob_samples <- rstan::extract(samples)$`strata_prob`
  mean_effect_samples <- rstan::extract(samples)$`mean_effect`
  treatment_effect_samples <- 0 * strata_prob_samples
  meta_data <- pstrata$PSobject$meta_data
  for (j in 1:ncol(treatment_effect_samples)) {
    g1 <- meta_data$SZDG[2 * j, 4] + 1
    g0 <- meta_data$SZDG[2 * j - 1, 4] + 1
    treatment_effect_samples[, j] <- 
      mean_effect_samples[, g1] - mean_effect_samples[, g0]
  }
  compare <- meta_data$strata
  strata_idx <- sapply(compare, function(x) which(x == meta_data$strata))
  contrast_effect_samples <- matrix(nrow = nrow(strata_prob_samples),
                                    ncol = length(strata_idx) * (length(strata_idx) - 1) / 2)
  cur_col_idx <- 1
  all_names <- c()
  for (i in 1:(length(strata_idx) - 1)) {
    for (j in (i + 1) : length(strata_idx)){
      all_names <- c(all_names,
                     paste0("Contrast: ", meta_data$strata[j], " over ", meta_data$strata[i], collapse = ' ')
      )
      contrast_effect_samples[, cur_col_idx] <- treatment_effect_samples[, strata_idx[j]] - treatment_effect_samples[, strata_idx[i]]
      cur_col_idx = cur_col_idx + 1
    }
  }
  contrast_effect_summary <- t(apply(contrast_effect_samples, 2, function(x) c(
    mean = mean(x), sd = sd(x),
    `2.5%` = unname(quantile(x, 0.025)), `25%` = unname(quantile(x, 0.25)),
    `median` = unname(quantile(x, 0.5)),
    `75%` = unname(quantile(x, 0.75)), `97.5%` = unname(quantile(x, 0.975))
  )))
  rownames(contrast_effect_summary) <- all_names
  return (contrast_effect_summary)
}


#' @method  summary PStrata_survival
#' @export summary.PStrata_survival
summary.PStrata_survival <- function(pstrata) {
  samples <- pstrata$post_samples
  strata_prob_samples <- rstan::extract(samples)$`strata_prob`
  mean_surv_prob_samples <- rstan::extract(samples)$`mean_surv_prob`
  treatment_effect_samples <- array(0, dim = c(dim(strata_prob_samples), dim(mean_surv_prob_samples)[3]))
  strata_prob_summary <- t(apply(strata_prob_samples, 2, function(x) c(
    mean = mean(x), sd = sd(x),
    `2.5%` = unname(quantile(x, 0.025)), `25%` = unname(quantile(x, 0.25)),
    `median` = unname(quantile(x, 0.5)),
    `75%` = unname(quantile(x, 0.75)), `97.5%` = unname(quantile(x, 0.975))
  )))
  mean_surv_prob <- apply(mean_surv_prob_samples, c(2,3), mean)
  lwr_surv_prob <- apply(mean_surv_prob_samples, c(2,3), quantile, 0.025)
  upr_surv_prob <- apply(mean_surv_prob_samples, c(2,3), quantile, 0.975)
  colnames(mean_surv_prob) <- round(psresult$PSobject$stan_data$time, 3)
  colnames(lwr_surv_prob) <- round(psresult$PSobject$stan_data$time, 3)
  colnames(upr_surv_prob) <- round(psresult$PSobject$stan_data$time, 3)
  meta_data <- pstrata$PSobject$meta_data
  rownames(strata_prob_summary) <- paste("S = ", meta_data$strata, sep = '')
  group_names <- rep(NA, nrow(mean_surv_prob))
  for (i in 1:nrow(meta_data$SZDG)){
    cur_g = meta_data$SZDG[i, 4] + 1
    cur_s = meta_data$SZDG[i, 1] + 1
    cur_z = meta_data$SZDG[i, 2]
    if (!is.na(group_names[cur_g])) {
      group_names[cur_g] <- 
        paste0("S = ", meta_data$strata[cur_s], ", Z = *")
    }
    else {
      group_names[cur_g] <- 
        paste0("S = ", meta_data$strata[cur_s], ", Z = ", cur_z)
    }
  }
  rownames(mean_surv_prob) <- group_names
  rownames(lwr_surv_prob) <- group_names
  rownames(upr_surv_prob) <- group_names
  
  for (j in 1:dim(treatment_effect_samples)[2]) {
    g1 <- meta_data$SZDG[2 * j, 4] + 1
    g0 <- meta_data$SZDG[2 * j - 1, 4] + 1
    treatment_effect_samples[, j, ] <- 
      mean_surv_prob_samples[, g1, ] - mean_surv_prob_samples[, g0, ]
  }
  mean_treatment_effect <- apply(treatment_effect_samples, c(2,3), mean)
  lwr_treatment_effect <- apply(treatment_effect_samples, c(2,3), quantile, 0.025)
  upr_treatment_effect <- apply(treatment_effect_samples, c(2,3), quantile, 0.975)
  rownames(mean_treatment_effect) <- paste("S = ", meta_data$strata, sep = '')
  rownames(upr_treatment_effect) <- paste("S = ", meta_data$strata, sep = '')
  rownames(lwr_treatment_effect) <- paste("S = ", meta_data$strata, sep = '')
  colnames(mean_treatment_effect) <- round(psresult$PSobject$stan_data$time, 3)
  colnames(upr_treatment_effect) <- round(psresult$PSobject$stan_data$time, 3)
  colnames(lwr_treatment_effect) <- round(psresult$PSobject$stan_data$time, 3)
  
  return (list(strata_prob = strata_prob_summary,
               survival_prob = list(
                 mean = mean_surv_prob,
                 lwr = lwr_surv_prob,
                 upr = upr_surv_prob
               ),
               treatment_effect = list(
                 mean = mean_treatment_effect,
                 lwr = lwr_treatment_effect,
                 upr = upr_treatment_effect
               )
          )
    )
}

#' @method  plot PStrata_survival
#' @export plot.PStrata_survival
plot.PStrata_survival <- function(pstrata, 
                         type = c("strata_prob", "surv_prob", "treatment_effect")) {
  samples <- pstrata$post_samples
  meta_data <- pstrata$PSobject$meta_data
  strata_prob_samples <- rstan::extract(samples)$`strata_prob`
  colnames(strata_prob_samples) <- meta_data$strata
  strata_prob_plot <- tidyr::pivot_longer(as.data.frame(strata_prob_samples), 
                                          everything(), names_to = 'Stratum',
                                          values_to = 'Probability')
  summary_res <- summary.PStrata_survival(pstrata)
  break_rows <- function(mat, label) {
    .names <- rownames(mat)
    lapply(.names, function(name) {
      data.frame(value = mat[name, ], time = pstrata$PSobject$stan_data$time, name = name, label = label)
    }) %>% {do.call(bind_rows, .)}
  }
  df_surv_prob <- bind_rows(
    break_rows(summary_res$survival_prob$mean, "mean"),
    break_rows(summary_res$survival_prob$upr, "upr"),
    break_rows(summary_res$survival_prob$lwr, "lwr")
  )
  
  df_treatment_effect <- bind_rows(
    break_rows(summary_res$treatment_effect$mean, "mean"),
    break_rows(summary_res$treatment_effect$upr, "upr"),
    break_rows(summary_res$treatment_effect$lwr, "lwr")
  )
  
  Gplot <- NULL
  if (type == "strata_prob" || type == 1) {
    Gplot <- ggplot2::ggplot(strata_prob_plot, ggplot2::aes(Probability)) +
      ggplot2::geom_density(ggplot2::aes(fill = Stratum), alpha = 0.4)
  }
  else if (type == "surv_prob" || type == 2) {
    Gplot <- ggplot2::ggplot(df_surv_prob %>% filter(label == "mean"), 
                             ggplot2::aes(time, value, color = name)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(aes(
        ymax = df_surv_prob %>% filter(label == "upr") %>% pull(value),
        ymin = df_surv_prob %>% filter(label == "lwr") %>% pull(value),
        fill = name
      ), alpha = 0.4) +
      ylab("Survial Probility") + xlab("Time") + labs(color = "Group", fill = "Group")
  }
  else if (type == "treatment_effect" || type == 3) {
    Gplot <- ggplot2::ggplot(df_treatment_effect %>% filter(label == "mean"), 
                             ggplot2::aes(time, value, color = name)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(aes(
        ymax = df_treatment_effect %>% filter(label == "upr") %>% pull(value),
        ymin = df_treatment_effect %>% filter(label == "lwr") %>% pull(value),
        fill = name
      ), alpha = 0.4) +
      ylab("Survial Probility") + xlab("Time") + labs(color = "Group", fill = "Group")
  }
  
  return (Gplot)
}
