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
