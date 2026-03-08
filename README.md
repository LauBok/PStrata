# PStrata

<!-- badges: start -->
<!-- badges: end -->

**PStrata** fits Bayesian principal stratification models using **Stan**. It supports a wide variety of models, priors, and assumptions to provide flexibility for causal inference in the presence of post-treatment confounding.

See Liu and Li (2023) [arXiv:2304.02740](https://arxiv.org/abs/2304.02740) for details.

## Installation

```r
# Install from GitHub
devtools::install_github("LauBok/PStrata")
```

## Overview

Let $Z$ denote the assigned treatment, $D$ denote the post-randomization intermediate variable, $X$ denote the covariates and $Y$ denote the outcome. Define $S = (D(0), D(1))$ as the principal stratum.

**PStrata** requires specification of two models:

- **S-model**: a multinomial model for the distribution of principal strata given covariates $$\log \frac{p(S = s\mid X)}{p(S = s_0 \mid X)} = X \beta^S.$$

- **Y-model**: a generalized linear model for the outcome specific to a principal stratum and treatment, given the covariates $$\mathbb{E}[Y \mid X, S, Z] = g^{-1}(X \beta_{s, z}^Y).$$

## Assumptions

- **Monotonicity**: assumes $D(0)\leq D(1)$, ruling out the existence of principal stratum $S = (1, 0)$. In **PStrata**, users explicitly specify all principal strata that are assumed to exist, allowing more flexibility than the standard monotonicity assumption.

- **Exclusion restriction (ER)**: for principal stratum $S$ where $D(0) = D(1)$, assumes $p(Y \mid X, S, Z) = p(Y \mid X, S)$.

## Example

### Normal outcome

Consider fitting a Bayesian principal stratification model on `sim_data_normal`, with intercept-only S-model and Y-model. The Y-model uses a Gaussian family with identity link. We assume monotonicity and exclusion restriction on both $S = (0, 0)$ and $S = (1, 1)$.

```r
library(PStrata)

obj <- PStrata(
    S.formula = Z + D ~ 1,
    Y.formula = Y ~ 1,
    Y.family = gaussian(link = "identity"),
    data = sim_data_normal,
    strata = c(n = "00*", c = "01", a = "11*"),
    prior_intercept = prior_normal(0, 1),
    prior_sigma = prior_inv_gamma(1),
    chains = 4, warmup = 500, iter = 1000
)
```

The `strata` argument specifies the assumed principal strata. The `*` suffix denotes that the exclusion restriction is applied to that stratum. Names (e.g., `n`, `c`, `a` for never-taker, complier, always-taker) are optional.

Print `obj` for an overview of estimated stratum proportions. Use `summary(obj)` for quantiles and confidence intervals.

```r
obj
summary(obj)
```

### Outcome estimation and contrasts

To obtain estimated mean effects by principal stratum and treatment arm:

```r
res <- PSOutcome(obj)
summary(res, "data.frame")
plot(res)
```

To compute stratum-specific treatment effects ($Y(1) - Y(0)$):

```r
cont <- PSContrast(res, Z = TRUE)
summary(cont, "data.frame")
plot(cont)
```

### Survival outcome

**PStrata** also supports survival outcomes with Cox proportional hazards models:

```r
obj <- PStrata(
    S.formula = Z + D ~ 1,
    Y.formula = survival(time, status, method = "Cox") ~ 1,
    data = sim_data_Cox,
    strata = c(n = "00*", c = "01", a = "11*"),
    prior_intercept = prior_normal(0, 1),
    chains = 4, warmup = 500, iter = 1000
)

res <- PSOutcome(obj)
cont <- PSContrast(res, Z = TRUE)
plot(cont)
```

## Key Functions

| Function | Description |
|---|---|
| `PStrata()` | Fit a principal stratification model |
| `PSOutcome()` | Extract estimated outcomes by stratum and treatment |
| `PSContrast()` | Compute contrasts (e.g., treatment effects) |
| `PSFormula()` | Create a formula object for PStrata |
| `PStrataInfo()` | Specify strata and assumptions |
| `PSObject()` | Build the full model specification |
| `make_stancode()` | Generate Stan code for inspection |
| `make_standata()` | Generate Stan data for inspection |

## Prior distributions

**PStrata** provides a set of prior distribution functions:

`prior_normal()`, `prior_t()`, `prior_cauchy()`, `prior_logistic()`, `prior_exponential()`, `prior_gamma()`, `prior_inv_gamma()`, `prior_chisq()`, `prior_inv_chisq()`, `prior_weibull()`, `prior_lasso()`, `prior_flat()`

## License

GPL (>= 2)
