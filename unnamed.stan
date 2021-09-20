functions {
    real inv_gaussian_lpdf(real x, real mu, real lambda) {
        real constant = log(lambda) / 2.0 - log(2 * pi()) / 2.0;
        real kernel = - 1.5 * log(x) - lambda * pow(x - mu, 2) / (2 * x * pow(mu, 2));
        return constant + kernel;
    }
    real survival_lpdf(real x, real mu, real theta1, real theta2, int censor) {
        real term1 = theta1 + theta2 + mu + (exp(theta1) - 1) * log(x);
        real term2 = exp(theta2 + mu) * pow(x, exp(theta1));
        return (1 - censor) * term1 - term2;
    }
}
data {
    int<lower=1> N;
    int<lower=0, upper=1> Z[N];
    int<lower=0, upper=1> D[N];
    real Y[N];
}
parameters {
    real par_1;
    real par_2;
    real<lower=0> par_3;
    real par_4;
    real<lower=0> par_5;
    real par_6;
    real<lower=0> par_7;
}
model {
    par_3 ~ inv_gamma(1, 1);
    par_5 ~ inv_gamma(1, 1);
    par_7 ~ inv_gamma(1, 1);
    for (n in 1:N) {
        real log_prob_0 = 0;
        real log_prob_1 = par_1 * 1;
        real log_prob_all = log_sum_exp({log_prob_0, log_prob_1});
        int length;
        if (Z[n] == 0 && D[n] == 0)
            length = 2;
        else if (Z[n] == 1 && D[n] == 0)
            length = 1;
        else if (Z[n] == 1 && D[n] == 1)
            length = 1;
        {
            real log_l[length];
            if (Z[n] == 0 && D[n] == 0) {
                // strata: 0, 1
                log_l[1] = log_prob_0 + normal_lpdf(Y[n] | par_2 * 1, par_3);
                log_l[2] = log_prob_1 + normal_lpdf(Y[n] | par_4 * 1, par_5);
            }
            else if (Z[n] == 1 && D[n] == 0) {
                // strata: 0
                log_l[1] = log_prob_0 + normal_lpdf(Y[n] | par_2 * 1, par_3);
            }
            else if (Z[n] == 1 && D[n] == 1) {
                // strata: 1
                log_l[1] = log_prob_1 + normal_lpdf(Y[n] | par_6 * 1, par_7);
            }
            target += log_sum_exp(log_l) - log_prob_all;
        }
    }
}
