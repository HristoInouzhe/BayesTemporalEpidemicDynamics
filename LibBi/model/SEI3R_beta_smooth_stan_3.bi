model SEI3R_beta_smooth { 
  const N = 2188017;      // Total population size
  const h = 1e-1;         // Integration time step for ODE solver
  const a_1 = 999993.424608, a_2 = 4.575392, a_3 = 1.0, a_4 = 1.0; // Dirichlet hyperparameters via Gamma
  const gamma = 0.1;                       // Recovery rate per stage
  const alpha_mean = 0.5, alpha_sd = 0.05; // Prior mean and std for incubation rate α

  noise n_transmission; // Noise term driving log-transmission dynamics
  
  state S, E, I1, I2, I3, R, Z;  // Compartments: Susceptible, Exposed, 3 Infectious stages, Recovered, Cumulative incidence
  state log_beta;                // Log of transmission multiplier
  state S0, E0, I0, R0;          // Initial values of compartments
  state total;

  obs Incidence;  // Observed incidence at each discrete time step

  param p_R0;     // Initial reproduction number multiplier
  param phi_inv;  // Dispersion parameter for observation model
  param sig;      // Volatility of log_beta diffusion process
  param alpha;    // Transition rate from E to I1 (1 / incubation period)
  param X1, X2, X3, X4; // Latent Gamma variables for initial conditions (used to create Dirichlet sample)

  sub parameter { 
    p_R0    ~ uniform(0.01, 1);            // Prior on reproduction number
    phi_inv ~ uniform(0.01, 1);         // Prior on dispersion (inverse scale)
    sig     ~ uniform(0.01, 1);           // Prior on noise intensity in log_beta
    alpha   ~ gaussian(mean = alpha_mean, std = alpha_sd); // Prior on incubation rate

    // Gamma sampling to induce a Dirichlet distribution over initial states
    X1 ~ gamma(a_1, 1);
    X2 ~ gamma(a_2, 1);
    X3 ~ gamma(a_3, 1);
    X4 ~ gamma(a_4, 1);

    total <- X1 + X2 + X3 + X4;

    // Normalize to population scale to define initial states
    S0 <- (X1 / total) * N;
    E0 <- (X2 / total) * N;
    I0 <- (X3 / total) * N;
    R0 <- (X4 / total) * N;
  } 

  sub initial { 
    S  <- S0;       // Initial susceptible
    E  <- E0;       // Initial exposed
    I1 <- I0;       // Initial infectious in stage 1
    I2 <- 0;        // Start with 0 in stage 2
    I3 <- 0;        // Start with 0 in stage 3
    R  <- R0;       // Initial recovered
    Z  <- E0 + I0 + R0; // Total non-susceptible population at t=0
  }   

  sub transition (delta = h) { // Transition dynamics with time step h
    n_transmission ~ wiener(); // Brownian motion for stochastic log_beta

    Z <- (t_now % 1 == 0 ? 0 : Z); // Reset cumulative incidence daily

    // ODE integration of SEI₃R model with stochastic transmission
    ode (alg = 'RK4(3)', h = h, atoler = 1e-2, rtoler = 1e-5) { 
      dS/dt = - (p_R0 * exp(log_beta)) * S * (I1 + I2 + I3) / N; 
      dE/dt = (p_R0 * exp(log_beta)) * S * (I1 + I2 + I3) / N - alpha * E;
      dI1/dt = alpha * E - 3 * gamma * I1;
      dI2/dt = 3 * gamma * I1 - 3 * gamma * I2; 
      dI3/dt = 3 * gamma * I2 - 3 * gamma * I3; 
      dR/dt = 3 * gamma * I3; 
      dZ/dt = (p_R0 * exp(log_beta)) * S * (I1 + I2 + I3) / N; // Daily incidence
      dlog_beta/dt = sig * n_transmission / h; // Diffusion for log_beta (scaled Brownian motion)
    } 
  } 

  sub observation { 
    Incidence ~ negbin(Z, 1 / phi_inv); // Negative binomial observation model
  } 
}

