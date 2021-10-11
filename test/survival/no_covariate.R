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
n <- nrow(data)

get_one <- function(log_p1, log_p2, log_p3, log_p4) {
  m <- max(log_p1, log_p2, log_p3, log_p4, na.rm = T)
  f <- function(x) if (is.na(x)) 0 else exp(x - m)
  return (which.max(rmultinom(1, 1, c(f(log_p1), f(log_p2), f(log_p3), f(log_p4)))))
}

sample_survival <- function(n, theta_1, exponential){
  # CDF: 1 - exp(-exp(theta_2)*t^(exp(theta_1)))
  u <- runif(n)
  t <- (-log(u) / exp(exponential))^exp(-theta_1)
  return (t)
}

## S
data$S <- sapply(
  1:n,
  function(i) get_one(
    log(0.3), 
    log(0.2),
    log(0.5),
    NA
  )
)
data$Z <- rbinom(n, 1, 0.5)
data$D <- ifelse(data$Z == 1, 
                 ifelse(data$S %in% c(2, 4), 1, 0), 
                 ifelse(data$S %in% c(3, 4), 1, 0))

data$Y <- ifelse(data$S == 1,
                 sample_survival(n, 1, 1), 
                 ifelse(data$S == 2,
                        sample_survival(n, 1, 4 - 1 * data$Z),
                        sample_survival(n, 1, 6))
)

data$C <- rbinom(n, 1, 0)

write.csv(data, "test/survival/data_no_covariate.csv", row.names = F)

PSObject(
  S.model = Z + D ~ 1,
  Y.model = Y + C ~ 1,
  Y.family = survival(),
  monotonicity = "default",
  ER = c('00')
) -> obj

result <- PStrata(
  S.formula = Z + D ~ 1,
  Y.formula = Y + C ~ 1,
  Y.family = survival(),
  data = data,
  monotonicity = "default",
  ER = c('00', '11'),
  prior_coefficient = prior_normal(0, 10),
  trunc = FALSE,
  chains = 1, warmup = 200, iter = 500
)

result
saveRDS(result, file = "test/survival/no_covariate_result.RDS")
