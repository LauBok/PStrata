# Manual (Version 0.0.1) Updated: 3/24/2022

The package **PStrata** fits Bayesian principal stratification models using **Stan**. A wide variety of models, priors and assumptions are allowed to provide flexibility.

------------------------------------------------------------------------

## Overview

Let $Z$ denote the assigned treatment, $D$ denote the post-randomization intermediate variable, $X$ denote the covariates and $Y$ denote the outcome. Define $S = (D(0), D(1))$ as the principal stratum.

**PStrata** requires specification of two models.

-   S-model: a multinomial model that specifies the distribution of principal strata given covariates $$\log \frac{p(S = s\mid X)}{p(S = s_0 \mid X)} = X \beta^S.$$

-   Y-model: a generalized linear model that specifies the distribution of outcome specific to a principal stratum and treatment, given the covariates $$\mathbb{E}[Y \mid X, S, Z] = g^{-1}(X \beta_{s, z}^Y).$$

--------------------------------------------------------------------------

## Assumption

- monotonicity: the monitonicity assumption assumes that $D(0)\leq D(1)$, and hence ruling out the existence of principal stratum $S = (1, 0)$. In **PStrata**, users explicitly specify all principal strata that are assumed to exist. This allows for more flexibility than the original monotonicity assumption.

- exclusion restriction (ER): the ER assumption for principal stratum $S$ assumes that if $D(0) = D(1)$, then $p(Y \mid X, S, Z) = p(Y \mid X, S)$.

--------------------------------------------------------------------------

## Example

Consider fitting a Bayesian principal stratification model on dataset `sim_data_normal`, with both S-model and Y-model involving only intercepts and no covariates. The Y-model is a linear model (gaussian distribution with link function being the identity function). Assume monotonicity and exclusion restriction on both $S = (0, 0)$ and $S = (1, 1)$.

We specify standard normal priors for the intercepts included in the S-model and Y-model. Also, the Y-model introduces a standard deviation parameter $\sigma$ on which we specify a standard inverse gamma prior.

The following script runs the Bayesian analysis (may take around 2 minutes).
```
obj <- PStrata(
    S.formula = Z + D ~ 1, # treatment + intermediate ~ covariates
    Y.formula = Y ~ 1,     # outcome ~ covariates
    Y.family = gaussian(link = "identity"),
    data = sim_data_normal,
    strata = c(n = "00*", c = "01", a = "11*"), # * for ER, names are optional
    prior_intercept = prior_normal(0, 1),
    prior_sigma = prior_inv_gamma(1),
    chains = 4, warmup = 500, iter = 1000 # optional parameters to pass to rstan()
)
```

Raw output of `obj` provides an overview of the estimated proportion for each strata and `summary(obj)` provides quantiles and confidence intervals.

To obtain the estimated mean effects of each principal stratum and treatment arm, run the following script.
```
res <- PSOutcome(obj)
summary(res, "data.frame") # returns a data.frame object. Other options: "array", "matrix"
plot(res)
```

To obtain contrasts, use function `PSContrast`. For example, to compare the stratum-specific treatment effects (defined by $Y(1) - Y(0)$), use the following script.
```
cont <- PSContrast(res, Z = TRUE)
summary(cont, "data.frame")
plot(cont)
```

