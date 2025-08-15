functions {
  // ---------------------------------------------
  // Linear B-spline interpolation
  // ---------------------------------------------
  vector bspline_interp(real t, int T_bs, vector ts_bs, int m, vector B_flat) {
    vector[m] bs;

    if (t <= ts_bs[1]) {
      for (i_bs in 1:m) bs[i_bs] = B_flat[i_bs];
      return bs;
    }
    if (t >= ts_bs[T_bs]) {
      for (j_bs in 1:m) bs[j_bs] = B_flat[(T_bs - 1) * m + j_bs];
      return bs;
    }
    int k_idx = 1;
    for (idx in 1:(T_bs - 1)) {
      if (t >= ts_bs[idx] && t < ts_bs[idx+1]) { k_idx = idx; break; }
    }
    real w = (t - ts_bs[k_idx]) / (ts_bs[k_idx+1] - ts_bs[k_idx]);
    for (l_bs in 1:m)
      bs[l_bs] = (1 - w) * B_flat[(k_idx - 1) * m + l_bs]
               + w       * B_flat[ k_idx      * m + l_bs];
    return bs;
  }

  // ---------------------------------------------
  // SIKR ODE: y = (S, I1..IK, R, C)
  // theta = [gamma, beta_1..beta_m]
  // x_i = [K, T_bs, m], x_r packs N, ts_bs, B_flat
  // ---------------------------------------------
  array[] real sikr_ode(real t,
                        array[] real y,
                        array[] real theta,
                        array[] real x_r,
                        array[] int  x_i) {
    int K    = x_i[1];
    int T_bs = x_i[2];
    int m    = x_i[3];

    real gamma = theta[1];
    vector[m] beta_vec;
    for (i in 1:m) beta_vec[i] = theta[i + 1];

    real N = x_r[1];
    vector[T_bs] ts_bs_local;
    for (i in 1:T_bs) ts_bs_local[i] = x_r[1 + i];
    int bs_length = T_bs * m;
    vector[bs_length] B_flat_local;
    for (i in 1:bs_length) B_flat_local[i] = x_r[1 + T_bs + i];

    vector[m] bs_vals = bspline_interp(t, T_bs, ts_bs_local, m, B_flat_local);
    real trans = exp(dot_product(beta_vec, bs_vals));

    int state_dim = K + 3;
    array[state_dim] real dydt;

    real S = y[1];
    real I_total = 0;
    for (i in 1:K) I_total += y[1 + i];   // I1..IK are y[2..K+1]

    // S'
    dydt[1] = - trans * S * I_total / N;

    // I1'
    dydt[2] = trans * S * I_total / N - K * gamma * y[2];

    // I2..IK'
    for (i in 2:K)
      dydt[1 + i] = K * gamma * y[i] - K * gamma * y[1 + i];

    // R'
    dydt[K + 2] = K * gamma * y[K + 1];

    // C' (incidence)
    dydt[K + 3] = trans * S * I_total / N;

    return dydt;
  }

  // numerically safe softplus (scale s)
  real softplus(real x, real s) {
    real z = x / s;
    if (z >  20) return x;
    if (z < -20) return s * exp(z);
    return s * log1p(exp(z));
  }
}

data {
  // observations
  int<lower=1> T_obs;
  real         t0;
  vector[T_obs] ts_obs;
  int           y[T_obs];
  vector[T_obs] eta;

  // model dims + spline
  int<lower=1> K;          // infectious compartments
  int<lower=1> m;          // # spline basis functions
  int<lower=1> T_bs;       // # spline knots
  vector[T_bs]  ts_bs;
  matrix[T_bs, m] B;
  matrix[m, m]  K_bs;

  // population & priors
  real<lower=0> N;
  vector<lower=0>[3] alpha_dir;   // Dirichlet prior for (S0, I0, R0)
  real prior_mean_gamma;
  real<lower=0> prior_sd_gamma;
  real<lower=0> prior_rate_phi_inv;

  // fixing switches
  int<lower=0,upper=1> fix_gamma;
  real<lower=0>       gamma_fixed;
  int<lower=0,upper=1> fix_phi_inv;
  real<lower=0>       phi_inv_fixed;

  int<lower=0,upper=m> K0;        // how many betas to fix
  vector[m]           beta_fixed;  // first K0 entries used

  // initial-state fixing (S0,I0 only; R0 is free unless you fix both)
  int<lower=0,upper=1> fix_S0;
  real<lower=0,upper=1> S0_fixed;
  int<lower=0,upper=1> fix_I0;
  real<lower=0,upper=1> I0_fixed;
}

transformed data {
  // pack for ODE
  int x_i[3];
  x_i[1] = K;
  x_i[2] = T_bs;
  x_i[3] = m;

  vector[1 + T_bs + T_bs * m] x_r;
  x_r[1] = N;
  {
    int pos = 1;
    for (r in 1:T_bs) { pos += 1; x_r[pos] = ts_bs[r]; }
    for (r in 1:T_bs) for (c in 1:m) { pos += 1; x_r[pos] = B[r, c]; }
  }

  // number of free initial-state slots among (S0, I0, R0)
  int sum_fix_theta = fix_S0 + fix_I0;  // R0 not directly fixable here
  int F_theta = 3 - sum_fix_theta;
}

parameters {
  // free simplex over the remaining initial-state slots (S0,I0,R0)
  simplex[F_theta] theta_free;

  // epidemiology & obs
  real<lower=0> gamma_free;
  real<lower=0> phi_inv_free;
  real<lower=0> tau;

  // free spline coefficients beyond first K0 fixed ones
  vector[m - K0] beta_free;
}

transformed parameters {
  // rebuild (S0, I0, R0) fractions
  vector[3] theta;
  {
    real rem_mass_raw = 1
                      - (fix_S0 ? S0_fixed : 0)
                      - (fix_I0 ? I0_fixed : 0);
    real rem_mass = rem_mass_raw < 0 ? 0 : rem_mass_raw;
    int idx = 1;
    for (j in 1:3) {
      if (j == 1 && fix_S0)       { theta[1] = S0_fixed; }
      else if (j == 2 && fix_I0)  { theta[2] = I0_fixed; }
      else { theta[j] = theta_free[idx] * rem_mass; idx += 1; }
    }
  }

  // effective params
  real gamma_eff   = fix_gamma   ? gamma_fixed   : gamma_free;
  real phi_inv_eff = fix_phi_inv ? phi_inv_fixed : phi_inv_free;
  real phi_eff     = 1 / phi_inv_eff;

  // full beta
  vector[m] beta_eff;
  if (K0 > 0) for (i in 1:K0) beta_eff[i] = beta_fixed[i];
  for (i in (K0 + 1):m) beta_eff[i] = beta_free[i - K0];

  // pack theta for ODE: [gamma, beta...]
  array[m + 1] real theta_array;
  theta_array[1] = gamma_eff;
  for (i in 1:m) theta_array[i + 1] = beta_eff[i];

  // initial conditions in counts
  real S0 = N * theta[1];
  real I0 = N * theta[2];
  real R0 = N * theta[3];

  array[K + 3] real y0;
  y0[1] = S0;
  y0[2] = I0;                       // I1
  for (i in 3:(K + 1)) y0[i] = 0;   // I2..IK
  y0[K + 2] = R0;
  y0[K + 3] = N - S0;               // cumulative = I0 + R0

  array[T_obs] real ts_obs_arr = to_array_1d(ts_obs);

  // solve ODE
  array[T_obs, K + 3] real ode_sol =
    integrate_ode_rk45(sikr_ode,
                       y0,
                       t0,
                       ts_obs_arr,
                       theta_array,
                       to_array_1d(x_r),
                       x_i);

  // predicted new cases between observations
  vector[T_obs] C_pred;
  C_pred[1] = ode_sol[1, K + 3] - y0[K + 3];
  for (t in 2:T_obs)
    C_pred[t] = ode_sol[t, K + 3] - ode_sol[t - 1, K + 3];

  // guard against tiny negatives from numerics
  for (t in 1:T_obs) C_pred[t] = softplus(C_pred[t], 5);

  // expose for model/genq
  // (reuse names to avoid extra copies)
}

model {
  // Dirichlet prior on (S0,I0,R0) free slots
  {
    if (F_theta > 1) {
      vector[F_theta] alpha_dir_free;
      int idx = 1;
      for (j in 1:3) {
        if ((j==1 && fix_S0) || (j==2 && fix_I0)) continue;
        alpha_dir_free[idx] = alpha_dir[j];
        idx += 1;
      }
      theta_free ~ dirichlet(alpha_dir_free);
    }
  }

  // priors
  if (!fix_gamma)   gamma_free   ~ normal(prior_mean_gamma, prior_sd_gamma);
  if (!fix_phi_inv) phi_inv_free ~ exponential(prior_rate_phi_inv);

  // spline roughness prior
  tau ~ inv_gamma(1, 0.005);
  target += -0.5 / tau * quad_form(K_bs, beta_eff); 

  // observation model
  for (t in 1:T_obs)
    y[t] ~ neg_binomial_2(eta[t] * C_pred[t], phi_eff);
}

generated quantities {
  vector[T_obs] y_rep;
  for (t in 1:T_obs)
    y_rep[t] = neg_binomial_2_rng(eta[t] * C_pred[t], phi_eff); 
}

