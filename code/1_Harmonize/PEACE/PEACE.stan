data {
    int<lower=1> I ; // number of data points
    int<lower=1> K; // number of features
    int<lower=2> B; // number of batches
    array[I,K] real y ; // observations
    array[I] int<lower=1, upper=B> batch; // batches
    array[I] real x; // age
}
transformed data {
    int Z = 2; // number of mixture components
    vector[K] m_y ;
    vector[K] s_y ;
    vector[B] batch_sd ;
    vector[B] batch_counter ;
    real mtplr = 1;
    vector[B] multiplier_batch_sd ;
    vector[K] multiplier_s_y ;
    vector[B-1] inv_batch_sd ;
    vector[K] inv_s_y ;
    real s_x ;
    vector[K] multiplier_s_y_s_x ;
    for (k in 1:K) {
        m_y[k] = mean(y[1:I,k]);
        s_y[k] = sd(y[1:I,k]);
    }
    for (b in 1:B) {
        batch_sd[b] = 0;
        batch_counter[b] = 0;
    }
    for (i in 1:I) {
        for (k in 1:K) {
            batch_sd[batch[i]] += square(y[i,k] - m_y[k]);
            batch_counter[batch[i]] += 1;
        }
    }
    batch_sd = sqrt(batch_sd ./ (batch_counter-1));
    multiplier_batch_sd = mtplr * batch_sd ;
    multiplier_s_y = mtplr * s_y ;
    inv_s_y = 1 ./ s_y ;
    inv_batch_sd = 1 ./ batch_sd[2:B];
    s_x = sd(x);
    multiplier_s_y_s_x = mtplr * s_y ./ s_x ;
}
parameters {
    array[B] real<lower=0, upper=1> phi;
    array[Z] vector[K] alpha;
    array[Z] vector[K] beta;
    array[Z] vector<lower=0>[K] sigma;
    array[B-1] vector[K] gamma;
    array[B-1] vector<lower=0>[K] delta;
    vector<lower=0>[B-1] tau;
}
transformed parameters {
    array[B-1] vector[K] delta_sigma1;
    array[B-1] vector[K] delta_sigma2;
    for (b in 1:(B-1)) {
        delta_sigma1[b] = delta[b] .* sigma[1];
        delta_sigma2[b] = delta[b] .* sigma[2];
    }
}
model {
    tau ~ exponential(inv_batch_sd);
    for (b in 1:(B-1)) {
        gamma[b] ~ normal(0, tau[b]);
        delta[b] ~ cauchy(0, 1);
    }
    for (z in 1:Z) {
        alpha[z] ~ normal(m_y, multiplier_s_y);
        beta[z] ~ normal(0, multiplier_s_y_s_x);
        sigma[z] ~ exponential(inv_s_y);
    }
    for (i in 1:I) {
        int this_batch = batch[i];
        if (this_batch==1) {
            target += log_mix(phi[this_batch],
                normal_lpdf(y[i] | alpha[1] + beta[1] * x[i], sigma[1]),
                normal_lpdf(y[i] | alpha[2] + beta[2] * x[i], sigma[2]));
        } else {
            target += log_mix(phi[this_batch],
                normal_lpdf(y[i] | alpha[1] + beta[1] * x[i] + gamma[this_batch-1], delta_sigma1[this_batch-1]),
                normal_lpdf(y[i] | alpha[2] + beta[2] * x[i] + gamma[this_batch-1], delta_sigma2[this_batch-1]));
        }
    }
}
generated quantities {
    array[I,K] real y_combat ;
    {
        vector[Z] r;
        array[B] real log_phi = log(phi);
        array[B] real log_1mphi = log1m(phi);
        for (i in 1:I) {
            vector[Z] tmp ;
            int this_batch = batch[i];
            vector[K] alpha1_beta_x ;
            vector[K] alpha2_beta_x ;
            if (this_batch==1) {
                y_combat[i] = y[i] ;
            } else {
                alpha1_beta_x = alpha[1] + beta[1] * x[i];
                alpha2_beta_x = alpha[2] + beta[2] * x[i];
                tmp[1] = log_phi[this_batch] +
                    normal_lpdf(y[i] | alpha1_beta_x + gamma[this_batch-1], delta_sigma1[this_batch-1]);
                tmp[2] = log_1mphi[this_batch] +
                    normal_lpdf(y[i] | alpha2_beta_x + gamma[this_batch-1], delta_sigma2[this_batch-1]);
                r = softmax(tmp) ;
                for (k in 1:K) {
                    y_combat[i,k] =
                        r[1] * ((y[i,k] - alpha1_beta_x[k] - gamma[this_batch-1,k]) /
                            delta[this_batch-1,k] + alpha1_beta_x[k]) +
                        r[2] * ((y[i,k] - alpha2_beta_x[k] - gamma[this_batch-1,k]) /
                            delta[this_batch-1,k] + alpha2_beta_x[k]);
                }
            }
        }
    }
}

