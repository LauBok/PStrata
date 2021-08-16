# Test 1: No covariates, 2 strata (n - 00, c - 01), ER for 00.
# July 6, 2021
# Models:
# P(S = 00) = 0.3
# P(S = 01) = 0.7
# Y | S = 00 ~ N(3, 1)
# Y | S = 01, Z ~ N(-1 - Z, 2)

set.seed(0)
df <- datasets::infert
N <- nrow(df)

# Generate Z, S, Y, D
df$Z <- rbinom(N, 1, 0.5)
df$S <- sample(c(0, 1), size = N, replace = T, prob = c(0.3, 0.7))
df$D <- ifelse(df$S == 3 | df$S + df$Z == 2, 1, 0)
df$Y[df$S == 0] <- rnorm(sum(df$S == 0), 3, 1)
df$Y[df$S == 1] <- rnorm(sum(df$S == 1), -1 - df$Z, 2)

# PS
result <- PStrata::PStrata(
  S.formula = Z + D ~ 1,
  Y.formula = Y ~ 1,
  Y.family = gaussian(),
  data = df,
  monotonicity = "strong",
  ER = "00",
  trunc = FALSE,
  chains = 1, warmup = 200, iter = 500
)

rstan::stan(paste0("unnamed", ".stan"), data = df, 
            chains = 1, warmup = 200, iter = 500
)

tryCatch(
  rstan::stan_model(file = "unnamed.stan"),
  error = function(e) 
  {
    print(e$message) # or whatever error handling code you want
  },
  message = function(e) 
  {
    print(e$message) # or whatever error handling code you want
  }
)

zz <- file("test1.txt", open = "wt")
sink(zz ,type = "output")
sink(zz, type = "message")
rstan::stan_model(file = "unnamed.stan")
#and to close connections
sink()
sink()

plot(result)

tmp_obj <- PSobject_survival(
  S.formula = Z + D ~ X1 + I(X2^2) + X3 * X4,
  Y.formula = Y ~ W1 + W2*W3,
  monotonicity = "strong",
  ER = "00",
  trunc = F
)
