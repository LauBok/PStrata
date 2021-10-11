#gitcreds::gitcreds_set()

# Test 1: No covariates, 2 strata (n - 00, c - 01), ER for 00.
# Aug 28, 2021
# Models:
# P(S = 00) = 0.3
# P(S = 01) = 0.7
# Y | S = 00 ~ N(3, 1)
# Y | S = 01, Z ~ N(-1 - Z, 0.5)

set.seed(0)

data <- read.csv("test/no_covariates/data1.csv")
n <- nrow(data)

get_one <- function(log_p1, log_p2, log_p3, log_p4) {
  m <- max(log_p1, log_p2, log_p3, log_p4, na.rm = T)
  f <- function(x) if (is.na(x)) 0 else exp(x - m)
  return (which.max(rmultinom(1, 1, c(f(log_p1), f(log_p2), f(log_p3), f(log_p4)))))
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
                     rnorm(n, 3, 1), 
                     rnorm(n, -1 - data$Z, 0.5))

write.csv(data, "test/no_covariates/data1.csv", row.names = F)


result <- PStrata::PStrata(
  S.formula = Z + D ~ 1,
  Y.formula = Y ~ 1,
  Y.family = gaussian(),
  data = data,
  monotonicity = "strong",
  ER = c('00'),
  trunc = FALSE,
  chains = 1, warmup = 200, iter = 500
)


plot(result)
