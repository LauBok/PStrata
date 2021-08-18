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
    int<lower=0, upper=1> C[N];
    real<lower=0> Y[N];
}
parameters {
    real par_1;
    real par_2;
    real par_3;
    real par_4;
    real par_5;
    real par_6;
    real par_7;
    real par_8;
    real par_9;
    real par_10;
}
model {
    par_3 ~ normal(0, 1);
    par_4 ~ normal(0, 1);
    par_6 ~ normal(0, 1);
    par_7 ~ normal(0, 1);
    par_9 ~ normal(0, 1);
    par_10 ~ normal(0, 1);
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
                log_l[1] = log_prob_0 + survival_lpdf(Y[n] | par_2 * 1, par_3, par_4, C[n]);
                log_l[2] = log_prob_1 + survival_lpdf(Y[n] | par_5 * 1, par_6, par_7, C[n]);
            }
            else if (Z[n] == 1 && D[n] == 0) {
                // strata: 0
                log_l[1] = log_prob_0 + survival_lpdf(Y[n] | par_2 * 1, par_3, par_4, C[n]);
            }
            else if (Z[n] == 1 && D[n] == 1) {
                // strata: 1
                log_l[1] = log_prob_1 + survival_lpdf(Y[n] | par_8 * 1, par_9, par_10, C[n]);
            }
            target += log_sum_exp(log_l) - log_prob_all;
        }
    }
}
