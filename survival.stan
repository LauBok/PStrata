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
    int<lower=0> DimS;
    int<lower=0> DimY;
    matrix[N, DimS] XS;
    matrix[N, DimY] XY;
    real<lower=0> Y[N];
}
parameters {
    real S_01_Intrcpt;
    vector[3] S_01_Coef;
    real S_11_Intrcpt;
    vector[3] S_11_Coef;
    real Y_00_0_Intrcpt;
    vector[2] Y_00_0_Coef;
    real Y_00_0_Theta;
    real Y_01_0_Intrcpt;
    vector[2] Y_01_0_Coef;
    real Y_01_0_Theta;
    real Y_01_1_Intrcpt;
    vector[2] Y_01_1_Coef;
    real Y_01_1_Theta;
    real Y_11_0_Intrcpt;
    vector[2] Y_11_0_Coef;
    real Y_11_0_Theta;
}
model {
    S_01_Coef ~ normal(0, 1);
    S_11_Coef ~ normal(0, 1);
    Y_00_0_Coef ~ normal(0, 1);
    Y_00_0_Theta ~ normal(0, 1);
    Y_01_0_Coef ~ normal(0, 1);
    Y_01_0_Theta ~ normal(0, 1);
    Y_01_1_Coef ~ normal(0, 1);
    Y_01_1_Theta ~ normal(0, 1);
    Y_11_0_Coef ~ normal(0, 1);
    Y_11_0_Theta ~ normal(0, 1);
    for (n in 1:N) {
        real log_prob_0 = 0;
        real log_prob_1 = S_01_Intrcpt + XS*S_01_Coef;
        real log_prob_3 = S_11_Intrcpt + XS*S_11_Coef;
        real log_prob_all = log_sum_exp({log_prob_0, log_prob_1, log_prob_3});
        int length;
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
                // strata: 0, 1
                log_l[1] = log_prob_0 + survival_lpdf(. | Y_00_0_Theta + Y_00_0_Intrcpt + XY*Y_00_0_Coef, Y_00_0_Theta, .);
                log_l[2] = log_prob_1 + survival_lpdf(. | Y_01_0_Theta + Y_01_0_Intrcpt + XY*Y_01_0_Coef, Y_01_0_Theta, .);
            }
            else if (Z[n] == 0 && D[n] == 1) {
                // strata: 3
                log_l[1] = log_prob_3 + survival_lpdf(. | Y_11_0_Theta + Y_11_0_Intrcpt + XY*Y_11_0_Coef, Y_11_0_Theta, .);
            }
            else if (Z[n] == 1 && D[n] == 0) {
                // strata: 0
                log_l[1] = log_prob_0 + survival_lpdf(. | Y_00_0_Theta + Y_00_0_Intrcpt + XY*Y_00_0_Coef, Y_00_0_Theta, .);
            }
            else if (Z[n] == 1 && D[n] == 1) {
                // strata: 1, 3
                log_l[1] = log_prob_1 + survival_lpdf(. | Y_01_1_Theta + Y_01_1_Intrcpt + XY*Y_01_1_Coef, Y_01_1_Theta, .);
                log_l[2] = log_prob_3 + survival_lpdf(. | Y_11_0_Theta + Y_11_0_Intrcpt + XY*Y_11_0_Coef, Y_11_0_Theta, .);
            }
            target += log_sum_exp(log_l) - log_prob_all;
        }
    }
}
