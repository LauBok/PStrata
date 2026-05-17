# PStrata

<!-- badges: start -->
<!-- badges: end -->

**PStrata** fits Bayesian principal stratification models using **Stan**. It
supports a wide variety of outcome families, link functions, priors, and
customizable monotonicity / exclusion-restriction assumptions for causal
inference in the presence of post-treatment confounding.

See Liu and Li (2023) [arXiv:2304.02740](https://arxiv.org/abs/2304.02740) for
methodological details.

## Installation

```r
# Install from GitHub
devtools::install_github("LauBok/PStrata")
```

**PStrata** compiles models with **Stan** via `rstan`, so a working C++
toolchain is required (Rtools on Windows, Xcode command-line tools on macOS).
See the
[RStan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
guide.

## Overview

Let `Z` denote the assigned treatment, `D` the post-randomization intermediate
variable, `X` the covariates and `Y` the outcome. The principal stratum is
`S = (D(0), D(1))`. **PStrata** specifies two models:

- **S-model** — a multinomial (softmax) model for the principal stratum given
  covariates: `log p(S=s | X) / p(S=s0 | X) = X beta^S`.
- **Y-model** — a generalized linear (or survival) model for the outcome,
  specific to a principal stratum and treatment arm:
  `E[Y | X, S, Z] = g^{-1}(X beta_{s,z}^Y)`.

### Assumptions

- **Monotonicity** — you explicitly list every principal stratum assumed to
  exist via the `strata` argument; any stratum you omit (e.g. defiers `"10"`)
  is ruled out. This is more flexible than a fixed monotonicity assumption.
- **Exclusion restriction (ER)** — for a stratum where `D(0) = D(1)`, the
  outcome distribution does not depend on `Z`. Strata under ER are named in
  the `ER` argument.

## Workflow

| Step | Function | Purpose |
|---|---|---|
| 1 | `PStrataModel()` | Specify the model symbolically (no data needed) |
| 2 | `fit()` | Build data, generate Stan code, run MCMC |
| 3 | `estimate()` | Extract posterior potential outcomes per stratum × arm |
| 4 | `contrast()` | Compute causal effects (e.g. `Y(1) - Y(0)`) |

`print()`, `summary()`, `plot()`, `diagnostics()`, and `stancode()` are
available on the fitted object.

## Example

### Normal outcome

Fit an intercept-only model on `sim_data_normal` (non-compliance: never-takers,
compliers, always-takers) under monotonicity, with exclusion restriction on
never-takers and always-takers.

```r
library(PStrata)

model <- PStrataModel(
  S.formula = Z + D ~ 1,
  Y.formula = Y ~ 1,
  Y.family  = gaussian(link = "identity"),
  strata    = c(n = "00", c = "01", a = "11"),
  ER        = c("n", "a"),
  prior_intercept = prior_normal(0, 1),
  prior_sigma     = prior_inv_gamma(1)
)
summary(model)   # algebraic description of the specified model

ps_fit <- fit(model, data = sim_data_normal,
              chains = 4, warmup = 500, iter = 1000)
ps_fit            # stratum proportions + mean outcome per group
summary(ps_fit)   # posterior intervals for all parameters
diagnostics(ps_fit)
```

Each stratum is a string of `D` values, one digit per treatment arm: `"00"` =
never-taker, `"01"` = complier, `"11"` = always-taker. Names (`n`, `c`, `a`)
are optional. Multiple post-randomization variables use `|` (e.g. `"00|01"`)
or the list-of-lists form.

### Potential outcomes and contrasts

```r
est <- estimate(ps_fit)        # E[Y(z) | S = s] by stratum and arm
summary(est, "data.frame")
plot(est)

ctr <- contrast(ps_fit, Z = TRUE)   # treatment effect Y(1) - Y(0) per stratum
summary(ctr, "data.frame")
plot(ctr)
```

### Survival outcome

Time-to-event outcomes use `Y.family = survival("Cox")` (Weibull-Cox) or
`survival("AFT")` (log-normal AFT). The outcome formula carries the event
indicator: `time + status ~ X`.

```r
model_s <- PStrataModel(
  S.formula = Z + D ~ 1,
  Y.formula = Y + delta ~ 1,
  Y.family  = survival("Cox"),
  strata    = c(n = "00", c = "01", a = "11"),
  ER        = c("n", "a"),
  prior_intercept = prior_normal(0, 1)
)
ps_fit_s <- fit(model_s, data = sim_data_Cox,
                chains = 4, warmup = 500, iter = 1000)

# Survival probability (or type = "RACE" for restricted average causal effect)
est_s <- estimate(ps_fit_s, type = "probability")
plot(est_s)

ctr_s <- contrast(ps_fit_s, Z = TRUE)
plot(ctr_s)
```

## Outcome families

`gaussian`, `binomial`, `Gamma`, `poisson`, `inverse.gaussian` (standard
`stats::family` objects with their usual links), plus `survival("Cox")` and
`survival("AFT")`.

> Binary-outcome principal stratification is weakly identified at small sample
> sizes / short chains; use enough data and MCMC iterations (e.g.
> `chains = 4`, `iter = 2000`) and check `diagnostics()` for convergence.

## Prior distributions

`prior_normal()`, `prior_t()`, `prior_cauchy()`, `prior_logistic()`,
`prior_lasso()`, `prior_exponential()`, `prior_gamma()`, `prior_inv_gamma()`,
`prior_chisq()`, `prior_inv_chisq()`, `prior_weibull()`, `prior_flat()`.

Assign them through the `PStrataModel()` arguments `prior_intercept`,
`prior_coefficient`, `prior_sigma`, `prior_alpha`, `prior_lambda`,
`prior_theta`.

## Key functions

| Function | Description |
|---|---|
| `PStrataModel()` | Specify a principal stratification model |
| `fit()` | Fit the model via Stan MCMC |
| `estimate()` | Posterior potential outcomes by stratum × arm |
| `contrast()` | Causal-effect contrasts (strata, arms, time) |
| `diagnostics()` | MCMC convergence diagnostics |
| `make_stancode()` / `stancode()` | Inspect the generated Stan program |
| `survival()` | Family constructor for time-to-event outcomes |

## Datasets

`sim_data_normal` (Gaussian outcome, non-compliance) and `sim_data_Cox`
(survival outcome) are bundled for illustration.

## License

GPL (>= 2)