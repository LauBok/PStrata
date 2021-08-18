# Test 1: No covariates, 2 strata (n - 00, c - 01), ER for 00.
# July 6, 2021
# Models:
# P(S = 00) = 0.3
# P(S = 01) = 0.7
# Y | S = 00 ~ N(3, 1)
# Y | S = 01, Z ~ N(-1 - Z, 2)

set.seed(0)
#df <- datasets::infert
#N <- nrow(df)
N <- 300
df <- data.frame(id = 1:N)

# Generate Z, S, Y, D
df$Z <- rbinom(N, 1, 0.5)
df$S <- sample(c(0, 1), size = N, replace = T, prob = c(0.3, 0.7))
df$D <- ifelse(df$S == 3 | df$S + df$Z == 2, 1, 0)
df$Y[df$S == 0] <- exp(rnorm(sum(df$S == 0), 3, 1))
df$Y[df$S == 1] <- exp(rnorm(sum(df$S == 1), -1 - 2 * df$Z[df$S == 1], .5))

# PS
result <- PStrata::PStrata(
  S.formula = Z + D ~ 1,
  Y.formula = Y + D ~ 1,
  Y.family = survival(),
  data = df,
  monotonicity = "strong",
  ER = c('00'),
  trunc = FALSE,
  chains = 1, warmup = 200, iter = 500
)


plot(result)
