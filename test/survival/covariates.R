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
                 sample_survival(n, 1, 0.5 + 0.02 * data$age), 
                 ifelse(data$S == 2,
                        sample_survival(n, 1, 2 - data$Z - 0.05 * data$age),
                        sample_survival(n, 1, 1 - 0.05 * data$age))
)

data$C <- rbinom(n, 1, 0)

write.csv(data, "test/survival/data_covariate1.csv", row.names = F)

result <- PStrata::PStrata(
  S.formula = Z + D ~ age,
  Y.formula = Y + C ~ age ,
  Y.family = survival(),
  data = data,
  monotonicity = "strong",
  ER = c('00'),
  prior_coefficient = prior_normal(0, 1),
  trunc = FALSE,
  chains = 1, warmup = 200, iter = 500
)

result