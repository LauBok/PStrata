data {
    int<lower=1> N;
    int<lower=0, upper=1> Z[N];
    int<lower=0, upper=1> D[N];
    real Y[N];
}
parameters {
    real par_1;
    real par_2;
    real par_3;
    real<lower=0> par_4;
    real par_5;
    real<lower=0> par_6;
    real par_7;
    real<lower=0> par_8;
}
model {
    par_1 ~ normal(0, 100);
    par_2 ~ normal(0, 100);
    par_3 ~ normal(0, 100);
    par_4 ~ inv_gamma(0.1, 0.1);
    par_5 ~ normal(0, 100);
    par_6 ~ inv_gamma(0.1, 0.1);
    par_7 ~ normal(0, 100);
    par_8 ~ inv_gamma(0.1, 0.1);
    for (n in 1:N) {
    real log_prob_0 = par_1;
    real log_prob_1 = par_2;
    real log_prob_all = log_sum_exp({log_prob_0, log_prob_1});
        int length;
        if (Z[n] == 0 && D[n] == 0)
            length = 2;
        if (Z[n] == 1 && D[n] == 0)
            length = 1;
        if (Z[n] == 1 && D[n] == 1)
            length = 1;
        {
            real log_l[length];
            if (Z[n] == 0 && D[n] == 0) {
                // strata: 0, 1
                log_l[1] = log_prob_0 + normal_lpdf(Y[n] | par_3, par_4);
                log_l[2] = log_prob_1 + normal_lpdf(Y[n] | par_5, par_6);
            }
            if (Z[n] == 1 && D[n] == 0) {
                // strata: 0
                log_l[1] = log_prob_0 + normal_lpdf(Y[n] | par_3, par_4);
            }
            if (Z[n] == 1 && D[n] == 1) {
                // strata: 1
                log_l[1] = log_prob_1 + normal_lpdf(Y[n] | par_7, par_8);
            }
            target += log_sum_exp(log_l) - log_prob_all;
        }
    }
}
