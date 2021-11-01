# Test 2: No covariates, 3 strata (n - 00, c - 01, a - 11), ER for 00 and 11.
# Aug 28, 2021
# Models:
# P(S = 00) = 0.3
# P(S = 01) = 0.2
# P(S = 11) = 0.5
# Y | S = 00 ~ survival(hazard = 2)
# Y | S = 01, Z ~ survival(hazard = 1.5 - 0.5 * Z)

set.seed(0)

data <- read.csv("test/data.csv")
data <- list()
#n <- nrow(data)
n <- 10000
get_one <- function(log_p1, log_p2, log_p3, log_p4) {
  m <- max(log_p1, log_p2, log_p3, log_p4, na.rm = T)
  f <- function(x) if (is.na(x)) 0 else exp(x - m)
  return (which.max(rmultinom(1, 1, c(f(log_p1), f(log_p2), f(log_p3), f(log_p4)))))
}

sample_survival <- function(n, theta, mu){
  # CDF: 1 - exp(-exp(theta_2)*t^(exp(theta_1)))
  u <- runif(n)
  t <- (-log(1 - u) / exp(mu))^exp(-theta)
  return (t)
}

plot_survival <- function(theta, mu, t_range = c(0, 1)) {
  t_points <- seq(t_range[1], t_range[2], length.out = 100)
  surv_prob <- exp(-exp(mu) * t_points ^ exp(theta))
  ggplot(data.frame(time = t_points, probability = surv_prob)) + 
    geom_line(aes(time, surv_prob))
}

## S
data$S <- sapply(
  1:n,
  function(i) get_one(
    log(0.3), 
    log(0.7),
    NA,
    NA
  )
)
data$Z <- rbinom(n, 1, 0.5)
data$D <- ifelse(data$Z == 1, 
                 ifelse(data$S %in% c(2, 4), 1, 0), 
                 ifelse(data$S %in% c(3, 4), 1, 0))

data$Y <- ifelse(data$S == 1,
                 sample_survival(n, 1, 0.3), 
                 sample_survival(n, 1, 2 - 0.6 * data$Z))

p1 <- plot_survival(1, 0.3)
p2 <- plot_survival(1, 1.4)
p3 <- plot_survival(1, 2)

data$C <- rbinom(n, 1, 0.2)

ggplot(as.data.frame(data)) + geom_density(aes(x = Y, color = as.factor(Z))) +
  facet_wrap(~as.factor(S))

write.csv(data, "test/survival/data_no_covariate.csv", row.names = F)

PSobject <- PSObject(
  S.formula = Z + D ~ 1,
  Y.formula = Y + C ~ 1,
  Y.family = survival(),
  data = as.data.frame(data),
  monotonicity = "strong",
  ER = c('00'),
  prior_intercept = prior_normal(0, 1),
  prior_coefficient = prior_normal(0, 1),
  trunc = F
)

PSsample <- PSSampling(PSobject, "wrong", chains = 1, warmup = 300, iter = 1000, refresh = 10)
PSsampleEx <- PSSampleEx(PSobject, PSsample)
PSsummary <- PSSummary.survival(PSsampleEx)
p <- plot(PSsummary, time = seq(0.01, 1, length.out = 20))
(p1 / p2 / p3) | p

result <- PStrata(
  S.formula = Z + D ~ 1,
  Y.formula = Y + C ~ 1,
  Y.family = survival(),
  data = data,
  monotonicity = "strong",
  ER = c('00'),
  prior_intercept = prior_normal(0, 1),
  prior_coefficient = prior_normal(0, 1),
  trunc = FALSE,
  chains = 1, warmup = 1000, iter = 3000
)

result
saveRDS(result, file = "test/survival/no_covariate_result.RDS")
