#gitcreds::gitcreds_set()

# Test 2: No covariates, 3 strata (n - 00, c - 01, a - 11), ER for 00 and 11.
# Aug 28, 2021
# Models:
# P(S = 00) = 0.3
# P(S = 01) = 0.2
# P(S = 11) = 0.5
# Y | S = 00 ~ N(3, 1)
# Y | S = 01, Z ~ N(-1 - Z, 0.5)
# Y | S = 11 ~ N(1, 2)

set.seed(0)

data <- read.csv("test/data.csv")
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
    log(0.2),
    NA,
    log(0.5)
  )
)
data$Z <- rbinom(n, 1, 0.5)
data$D <- ifelse(data$Z == 1, 
                     ifelse(data$S %in% c(2, 4), 1, 0), 
                     ifelse(data$S %in% c(3, 4), 1, 0))
data$Y <- ifelse(data$S == 1,
                     rnorm(n, 3, 1), 
                 ifelse(data$S == 2,
                     rnorm(n, -1 - data$Z, 0.5),
                     rnorm(n, 1, 2))
)

write.csv(data, "test/no_covariates/data2.csv", row.names = F)


result <- PStrata(
  S.formula = Z + D ~ 1,
  Y.formula = Y ~ 1,
  Y.family = gaussian(),
  data = data,
  strata = c(0, 1, 3),
  ER = c(0, 3),
  chains = 2, warmup = 200, iter = 500
)

result
plot(result)


# if misspecified ER

result2 <- PStrata(
  S.formula = Z + D ~ 1,
  Y.formula = Y ~ 1,
  Y.family = gaussian(),
  data = data,
  monotonicity = "default",
  ER = c('00'),
  trunc = FALSE,
  chains = 1, warmup = 500, iter = 2000
)

plot(result2)
