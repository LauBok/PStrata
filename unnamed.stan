functions {
    real inv_gaussian_lpdf(real x, real mu, real lambda) {
        real constant = log(lambda) / 2.0 - log(2 * pi()) / 2.0;
        real kernel = - 1.5 * log(x) - lambda * pow(x - mu, 2) / (2 * x * pow(mu, 2));
        return constant + kernel;
    }
    real survival_lpdf(real x, real mu, real theta, int censor) {
        real term1 = theta + mu + (exp(theta) - 1) * log(x);
        real term2 = exp(mu) * pow(x, exp(theta));
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
    real S_01_Intrcpt;
    real Y_00_0_Intrcpt;
    real Y_00_0_Theta;
    real Y_01_0_Intrcpt;
    real Y_01_0_Theta;
    real Y_01_1_Intrcpt;
    real Y_01_1_Theta;
}
model {
    S_01_Intrcpt ~ normal(0, 1);
    Y_00_0_Intrcpt ~ normal(0, 1);
    Y_00_0_Theta ~ normal(0, 1);
    Y_01_0_Intrcpt ~ normal(0, 1);
    Y_01_0_Theta ~ normal(0, 1);
    Y_01_1_Intrcpt ~ normal(0, 1);
    Y_01_1_Theta ~ normal(0, 1);
    for (n in 1:N) {
        real log_prob_0 = 0;
        real log_prob_1 = S_01_Intrcpt;
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
                log_l[1] = log_prob_0 + survival_lpdf(Y[n] | Y_00_0_Theta + Y_00_0_Intrcpt, Y_00_0_Theta, C[n]);
                log_l[2] = log_prob_1 + survival_lpdf(Y[n] | Y_01_0_Theta + Y_01_0_Intrcpt, Y_01_0_Theta, C[n]);
            }
            else if (Z[n] == 1 && D[n] == 0) {
                // strata: 0
                log_l[1] = log_prob_0 + survival_lpdf(Y[n] | Y_00_0_Theta + Y_00_0_Intrcpt, Y_00_0_Theta, C[n]);
            }
            else if (Z[n] == 1 && D[n] == 1) {
                // strata: 1
                log_l[1] = log_prob_1 + survival_lpdf(Y[n] | Y_01_1_Theta + Y_01_1_Intrcpt, Y_01_1_Theta, C[n]);
            }
            target += log_sum_exp(log_l) - log_prob_all;
        }
    }
}
