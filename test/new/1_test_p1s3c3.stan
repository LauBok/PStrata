
data {
    int<lower=0> N;
    int<lower=0> PS;
    int<lower=0> PG;
    int<lower=0, upper=1> Z[N];
    int<lower=0, upper=1> D[N];
    int<lower=0, upper=1> Y[N];
    matrix[N, PS] XS;
    matrix[N, PG] XG;
}

transformed data {
    int S[4];
    S[1] = 1;
    S[2] = 2;
    S[3] = 2;
    S[4] = 3;
}

parameters {
    matrix[2, PS] beta_S;
    matrix[4, PG] beta_G;
}

model {
    beta_S[:, 1] ~ normal(0, 1);
    beta_G[:, 1] ~ normal(0, 1);
    to_vector(beta_S[:, 2:PS]) ~ normal(0.000000, 1.000000);
    to_vector(beta_G[:, 2:PG]) ~ normal(0.000000, 1.000000);
    for (n in 1:N) {
        int length;
        real log_prob[3];
        log_prob[1] = 0;
        for (s in 2:3) {
            log_prob[s] = XS[n] * beta_S[s - 1]';
        }
        if (Z[n] == 0 && D[n] == 0)
            length = 2;
        else if (Z[n] == 0 && D[n] == 1)
            length = 1;
        else if (Z[n] == 1 && D[n] == 0)
            length = 1;
        else if (Z[n] == 1 && D[n] == 1)
            length = 2;
        {
            real log_l[length];
            if (Z[n] == 0 && D[n] == 0) {
                // strata: 0 1
                log_l[1] = log_prob[1] + bernoulli_lpmf(Y[n] | inv_logit(XG[n] * beta_G[1]'));
                log_l[2] = log_prob[2] + bernoulli_lpmf(Y[n] | inv_logit(XG[n] * beta_G[2]'));
            }
            else if (Z[n] == 0 && D[n] == 1) {
                // strata: 2
                log_l[1] = log_prob[3] + bernoulli_lpmf(Y[n] | inv_logit(XG[n] * beta_G[4]'));
            }
            else if (Z[n] == 1 && D[n] == 0) {
                // strata: 0
                log_l[1] = log_prob[1] + bernoulli_lpmf(Y[n] | inv_logit(XG[n] * beta_G[1]'));
            }
            else if (Z[n] == 1 && D[n] == 1) {
                // strata: 1 2
                log_l[1] = log_prob[2] + bernoulli_lpmf(Y[n] | inv_logit(XG[n] * beta_G[3]'));
                log_l[2] = log_prob[3] + bernoulli_lpmf(Y[n] | inv_logit(XG[n] * beta_G[4]'));
            }
            target += log_sum_exp(log_l) - log_sum_exp(log_prob);
        }
    }
}

generated quantities {
    vector[4] mean_effect;
    {
        matrix[N, 4] expected_mean = XG * beta_G';
        matrix[N, 3] log_prob;
        vector[3] denom;
        vector[4] numer;
        log_prob[:, 1] = rep_vector(0, N);
        log_prob[:, 2:3] = XS * beta_S';
        for (n in 1:N) {
            log_prob[n] -= log_sum_exp(log_prob[n]);
        }
        for (s in 1:3) denom[s] = mean(exp(log_prob[:, s]));
        for (g in 1:4) {
            numer[g] = mean(expected_mean[:, g] .* exp(log_prob[:, S[g]]));
            mean_effect[g] = numer[g] / denom[S[g]];
        }
    }
}

