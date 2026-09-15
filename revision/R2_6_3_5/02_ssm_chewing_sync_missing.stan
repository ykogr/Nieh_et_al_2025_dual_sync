data {
  int<lower=1> N;           // Total Data Points（obs+mis）
  int<lower=0> N_obs;       // Number of observations
  int<lower=0> N_mis;       // Number of missing values
  int<lower=1> N_pair;
  int<lower=1> N_session;
  int<lower=1> N_bin;

  int<lower=1,upper=N> obs_idx[N_obs];
  int<lower=1,upper=N> mis_idx[N_mis];
  int<lower=1,upper=N_pair> pair_id[N];
  int<lower=1,upper=N_session> session_id[N];
  int<lower=1,upper=N_bin> bin[N];
  int<lower=0,upper=1> is_real[N];  // real=1, shuffled=0

  vector[N_obs] y_obs;
}
parameters {
  vector[N_mis] y_mis;    // Missing values

  real mu0;
  vector[N_pair] mu_pair;
  real<lower=0> sigma_pair;

  matrix[N_pair, N_session] mu_sess_init;
  real<lower=0> sigma_sess_init;

  vector[N_bin] mu_state;
  real<lower=0> sigma_mu;

  real beta0;
  vector<lower=-pi()/2, upper=pi()/2>[N_bin-1] beta_tan;
  real<lower=0> sigma_beta;

  real<lower=0> sigma_obs;
}
transformed parameters {
  vector[N_bin] beta_state;
  beta_state[1] = beta0;
  for (b in 2:N_bin)
    beta_state[b] = beta_state[b-1] + sigma_beta * tan(beta_tan[b-1]);
}
model {
  vector[N] y_full;
  {
    int o = 1;
    int m = 1;
    for (i in 1:N) {
      if (o <= N_obs && i == obs_idx[o]) {
        y_full[i] = y_obs[o];
        o += 1;
      } else {
        y_full[i] = y_mis[m];
        m += 1;
      }
    }
  }
  mu0 ~ normal(0, 1.5);
  sigma_pair ~ normal(0, 1)T[0,];
  sigma_sess_init ~ normal(0, 1)T[0,];
  sigma_mu ~ normal(0, 1)T[0,];
  sigma_beta ~ normal(0, 1)T[0,];
  sigma_obs ~ normal(0, 1)T[0,];

  mu_pair ~ normal(mu0, sigma_pair);
  to_vector(mu_sess_init) ~ normal(0, sigma_sess_init);

  mu_state[1] ~ normal(0, sigma_mu);
  for (b in 2:N_bin)
    mu_state[b] ~ normal(mu_state[b-1], sigma_mu);

  beta0 ~ normal(0, 1);
  beta_tan ~ uniform(-pi()/2, pi()/2);  // robust RW difference

  for (n in 1:N) {
    real mu = mu_pair[pair_id[n]]
            + mu_sess_init[pair_id[n], session_id[n]]
            + mu_state[bin[n]]
            + beta_state[bin[n]] * is_real[n];
    y_full[n] ~ normal(mu, sigma_obs);
  }
}
generated quantities {
  vector[N] y_full_out;
  {
    int o = 1;
    int m = 1;
    for (i in 1:N) {
      if (o <= N_obs && i == obs_idx[o]) {
        y_full_out[i] = y_obs[o];
        o += 1;
      } else {
        y_full_out[i] = y_mis[m];
        m += 1;
      }
    }
  }
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    real mu = mu_pair[pair_id[n]]
            + mu_sess_init[pair_id[n], session_id[n]]
            + mu_state[bin[n]]
            + beta_state[bin[n]] * is_real[n];
    y_rep[n] = normal_rng(mu, sigma_obs);
    log_lik[n] = normal_lpdf(y_full_out[n] | mu, sigma_obs);
  }
}
