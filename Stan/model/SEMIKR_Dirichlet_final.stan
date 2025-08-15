functions {
  // ------------------------------------------------------------
  // Linear B-spline interpolation: given time t, knots ts_bs,
  // basis evaluations B_flat (flattened T_bs×m matrix) and
  // number of basis functions m, return the m-vector of basis
  // weights at time t.
  // ------------------------------------------------------------
  vector bspline_interp(real t, int T_bs, vector ts_bs, int m, vector B_flat) {
    vector[m] bs;  // output basis values at t

    // if t is before the first knot, use the first row of B_flat
    if (t <= ts_bs[1]) {
      for (i_bs in 1:m)
        bs[i_bs] = B_flat[i_bs];
      return bs;
    }
    // if t is after the last knot, use the last row of B_flat
    if (t >= ts_bs[T_bs]) {
      for (j_bs in 1:m)
        bs[j_bs] = B_flat[(T_bs - 1) * m + j_bs];
      return bs;
    }
    // otherwise find the interval [ts_bs[k_idx], ts_bs[k_idx+1]) containing t
    int k_idx = 1;
    for (idx in 1:(T_bs - 1)) {
      if (t >= ts_bs[idx] && t < ts_bs[idx+1]) {
        k_idx = idx;
        break;
      }
    }
    // compute linear interpolation weight w ∈ [0,1]
    real w = (t - ts_bs[k_idx]) / (ts_bs[k_idx+1] - ts_bs[k_idx]);
    // interpolate between the two rows of B_flat
    for (l_bs in 1:m)
      bs[l_bs] = (1 - w) * B_flat[(k_idx - 1) * m + l_bs]
                + w       * B_flat[ k_idx      * m + l_bs];
    return bs;
  }

  // ------------------------------------------------------------
  // SEMIKR ODE system: state y, parameters theta = [alpha,gamma,β₁…βₘ],
  // data x_r,x_i.  Returns dy/dt.
  // ------------------------------------------------------------
  array[] real semikr_ode(real t,
                          array[] real y,
                          array[] real theta,
                          array[] real x_r,
                          array[] int  x_i) {
    // unpack dimensions
    int M    = x_i[1];   // # exposed compartments
    int K    = x_i[2];   // # infected compartments
    int T_bs = x_i[3];   // # spline knots
    int m    = x_i[4];   // # spline basis functions

    // extract alpha, gamma
    real alpha = theta[1];
    real gamma = theta[2];
    // extract β vector of length m
    vector[m] beta_vec;
    for (i in 1:m)
      beta_vec[i] = theta[i + 2];

    // reconstruct N, ts_bs, B_flat from x_r
    real N = x_r[1];  
    vector[T_bs] ts_bs_local;
    for (i in 1:T_bs)
      ts_bs_local[i] = x_r[1 + i];
    int bs_length = T_bs * m;
    vector[bs_length] B_flat_local;
    for (i in 1:bs_length)
      B_flat_local[i] = x_r[1 + T_bs + i];

    // evaluate time-varying transmission: exp(β·B(t))
    vector[m] bs_vals = bspline_interp(t, T_bs, ts_bs_local, m, B_flat_local);
    real trans = exp(dot_product(beta_vec, bs_vals));

    // dimension of the ODE state: S + M·E + K·I + R + cumulative cases
    int state_dim = M + K + 3;
    array[state_dim] real dydt;

    // compute total infected I_total = sum of the K infected compartments
    real S_val  = y[1];
    real I_total = 0;
    for (i in 1:K)
      I_total += y[M + i + 1];

    // SE...I...R transitions:
    dydt[1] = - trans * S_val * I_total / N;  // S′ = −β S I / N
    // ---------- E-stage (only if M>0) ----------
    real inflow_I1;                       // flow into the first I-compartment

    if (M > 0) {                    // usual SEMIKR dynamics
        dydt[2] = trans * S_val * I_total / N  // E₁′
                    - M * alpha * y[2];

        for (i in 2:M)                // E₂ … E_M
            dydt[i+1] = M * alpha * y[i] - M * alpha * y[i+1];

        inflow_I1 = M * alpha * y[M+1];
    } else {                         // SIKR: infections jump straight to I₁
        inflow_I1 = trans * S_val * I_total / N;
    }
    // I₁ equation:
    dydt[M+2] = inflow_I1 - K * gamma * y[M+2];
    // I₂...I_K exchange:
    for (i in 2:K)
      dydt[M + i + 1] = K * gamma * y[M + i] - K * gamma * y[M + i + 1];
    // R′ and cumulative cases:
    dydt[M + K + 2] = K * gamma * y[M + K + 1];             // R′
    dydt[M + K + 3] = trans * S_val * I_total / N;          // new cumulative cases

    return dydt;
  }

  real softplus(real x, real s) {
    real z = x / s;
    if (z >  20) return x;                     // for large x, softplus ≈ x
    if (z < -20) return s * exp(z);            // for very negative, use exp safely
    return s * log1p(exp(z));                  // otherwise the usual formula
  }
}

data {
  // ─────────────── Observation data ───────────────
  int<lower=1> T_obs;      // # observed time points
  real         t0;         // initial time
  vector[T_obs] ts_obs;    // observation times
  int           y[T_obs];  // observed counts
  vector[T_obs] eta;       // under report level

  // ─────── Model dimensions and spline data ───────
  int<lower=0> M;          // E compartments (0 ⇒ no exposed stage)
  int<lower=1> K;          // I compartments
  int<lower=1> m;          // number of spline basis functions
  int<lower=1> T_bs;       // # spline knots
  vector[T_bs]  ts_bs;     // spline knot times
  matrix[T_bs, m] B;       // spline basis evals
  matrix[m, m]  K_bs;      // spline roughness penalty

  // ─────────── Epidemic & prior constants ───────────
  real<lower=0> N;                 // population
  vector<lower=0>[4] alpha_dir;    // Dirichlet hyperparams for θ = (S₀,E₀,I₀,R₀)
  real prior_mean_alpha;           // prior for alpha
  real<lower=0> prior_sd_alpha;
  real prior_mean_gamma;           // prior for gamma
  real<lower=0> prior_sd_gamma;
  real<lower=0> prior_rate_phi_inv; // prior for overdispersion

  // ────────── Fixing switches (transmission) ──────────
  int<lower=0,upper=1> fix_alpha;
  real<lower=0>       alpha_fixed;
  int<lower=0,upper=1> fix_gamma;
  real<lower=0>       gamma_fixed;
  int<lower=0,upper=1> fix_phi_inv;
  real<lower=0>       phi_inv_fixed;
  int<lower=0,upper=m> K0;         // how many β’s to fix
  vector[m]         beta_fixed;   // those β values

  // ────────── Fixing switches (initial states) ──────────
  int<lower=0,upper=1> fix_S0;
  real<lower=0,upper=1>   S0_fixed;  // fraction of N for S₀
  int<lower=0,upper=1> fix_E0;
  real<lower=0,upper=1>   E0_fixed;  // fraction of N for E₀
  int<lower=0,upper=1> fix_I0;
  real<lower=0,upper=1>   I0_fixed;  // fraction of N for I₀
}

transformed data {
  // ─── pack ODE data for integrate_ode_rk45 ───
  // Declare x_i as a data-only integer array.
  int x_i[4];
  x_i[1] = M;
  x_i[2] = K;
  x_i[3] = T_bs;
  x_i[4] = m;
    
  // Build x_r as a data-only vector from N, ts_bs, and B.
  vector[1 + T_bs + T_bs * m] x_r;
  x_r[1] = N;
  {
    int pos = 1;
    for (dummy3 in 1:T_bs) {
      pos += 1;
      x_r[pos] = ts_bs[dummy3];
    }
    for (dummy4 in 1:T_bs) {
      for (dummy5 in 1:m) {
        pos += 1;
        x_r[pos] = B[dummy4, dummy5];
      }
    }
  }

  // ─── compute how many θ‐slots remain free ───
  int sum_fix_theta = fix_S0 + fix_E0 + fix_I0;
  int F_theta       = 4 - sum_fix_theta;  // number of free proportions
}

parameters {
  // ─── free simplex of remaining initial proportions ───
  simplex[F_theta] theta_free;

  // ─── free versions of α, γ, φ⁻¹, spline‐smoothness τ ───
  real<lower=0> alpha_free;
  real<lower=0> gamma_free;
  real<lower=0> phi_inv_free;
  real<lower=0> tau;

  // ─── free spline coefficients beyond the first K0 ───
  vector[m - K0] beta_free;
}

transformed parameters {
  // ─── rebuild full θ = (S₀,E₀,I₀,R₀) as a length‐4 vector ───
  vector[4] theta;
  {
    // rem_mass = leftover fraction after fixed slots
    real rem_mass_raw = 1
      - (fix_S0 ? S0_fixed : 0)
      - (fix_E0 ? E0_fixed : 0)
      - (fix_I0 ? I0_fixed : 0);
    real rem_mass = rem_mass_raw < 0 ? 0 : rem_mass_raw;  // ensure non‐negativity
    int idx = 1;
    for (j in 1:4) {
      if (j == 1 && fix_S0 == 1) {
        theta[1] = S0_fixed;  // S₀ locked in
      } else if (j == 2 && fix_E0 == 1) {
        theta[2] = E0_fixed;  // E₀ locked in
      } else if (j == 3 && fix_I0 == 1) {
        theta[3] = I0_fixed;  // I₀ locked in
      } else {
        // free slot: take share = theta_free[idx] × remaining mass
        theta[j] = theta_free[idx] * rem_mass;
        idx += 1;
      }
    }
  }

  // ─── effective transmission & overdispersion parameters ───
  real alpha_eff   = fix_alpha   ? alpha_fixed   : alpha_free;
  real gamma_eff   = fix_gamma   ? gamma_fixed   : gamma_free;
  real phi_inv_eff = fix_phi_inv ? phi_inv_fixed : phi_inv_free;
  real phi_eff     = 1 / phi_inv_eff;  // φ = 1/(φ⁻¹)

  // ─── rebuild full β of length m ───
  vector[m] beta_eff;
  if (K0 > 0) {
    // when K0==0 this entire block is skipped → no 1:0 loop
    for (i in 1:K0) {
      beta_eff[i] = beta_fixed[i];
    }
  }

  // 2) fill the remaining β’s from the free vector
  for (i in (K0 + 1):m) {
    // when K0==m this loop becomes (m+1):m, which also runs 0 times
    beta_eff[i] = beta_free[i - K0];
  }

  // ─── pack for ODE solver ───
  array[m+2] real theta_array;
  theta_array[1] = alpha_eff;
  theta_array[2] = gamma_eff;
  for (i in 1:m)
    theta_array[i+2] = beta_eff[i];

  // ─── initial compartments (in counts) ───
  real S0 = N * theta[1];  // S₀ = fraction × N
  real E0 = N * theta[2];
  real I0 = N * theta[3];
  // R₀ = leftover individuals
  real R0 = N * theta[4];

  // ─── set up and solve ODE ───
  array[M+K+3] real y0;
  y0[1] = S0;
  if (M > 0) {              // only when an E stage exists
    y0[2] = E0;
    for (i in 3:(M+1))
      y0[i] = 0;            // E₂ … E_M
  }
  y0[M+2] = I0;
  for (i in (M+3):(M+K+1)) y0[i] = 0;      // I₂…I_K
  y0[M+K+2] = R0;
  y0[M+K+3] = N - S0;                    // cumulative cases

  array[M+K+3] real y0_arr = to_array_1d(y0);
  array[T_obs]   real ts_obs_arr = to_array_1d(ts_obs);

  array[T_obs, M+K+3] real ode_sol =
    integrate_ode_rk45(semikr_ode,
                       y0_arr,
                       t0,
                       ts_obs_arr,
                       theta_array,
                       to_array_1d(x_r),
                       x_i);

  // ─── derive new‐case predictions C_pred[t] ───
  vector[T_obs] C_pred;
  C_pred[1] = ode_sol[1, M+K+3] - y0[M+K+3];
  for (t in 2:T_obs)
    C_pred[t] = ode_sol[t, M+K+3] - ode_sol[t-1, M+K+3];
  for (t in 1:T_obs)
    // C_pred[t] = fmax(C_pred[t], 1e-6);  // ensure non‐negativity of C_pred
    C_pred[t] = softplus(C_pred[t], 5);  // ensure non‐negativity of C_pred
}

model {
  // ─── Dirichlet prior on the free θ components ───
  {  
    if (F_theta > 1) {
      vector[F_theta] alpha_dir_free;
      int idx = 1;
      for (j in 1:4) {
         if ((j==1 && fix_S0) || (j==2 && fix_E0) || (j==3 && fix_I0))
           continue;
         alpha_dir_free[idx] = alpha_dir[j];
         idx += 1;
      }
      theta_free ~ dirichlet(alpha_dir_free);
    }
  }

  // ─── priors on transmission & overdispersion ───
  if (fix_alpha == 0)
    alpha_free     ~ normal(prior_mean_alpha, prior_sd_alpha);
  if (fix_gamma == 0)
    gamma_free     ~ normal(prior_mean_gamma, prior_sd_gamma);
  if (fix_phi_inv == 0)
    phi_inv_free   ~ exponential(prior_rate_phi_inv);

  // ─── roughness prior on spline coefficients ───
  tau             ~ inv_gamma(1, 0.005);
  target += -0.5 / tau * quad_form(K_bs, beta_eff);

  // ─── observation model ───
  for (t in 1:T_obs)
    y[t] ~ neg_binomial_2(eta[t] * C_pred[t], phi_eff);
}

generated quantities {
  vector[T_obs] y_rep;
  // posterior predictive cases
  for (t in 1:T_obs)
    y_rep[t] = neg_binomial_2_rng(eta[t] * C_pred[t], phi_eff);
}

