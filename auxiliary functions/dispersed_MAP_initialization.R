# Dispersed initialisation around the MAP for Stan models
# ------------------------------------------------------
# This helper creates multiple well‑spread starting points for an
# rstan::sampling() run.  The workflow is:
#   1. Draw n_inits_0 seeds from the *priors*.
#   2. Run a short LBFGS from each seed → approximate MAP.
#   3. Pick the best MAP and *jitter* it to obtain n_inits proposals.
#   4. Score every proposal with the exact log‑posterior.
#   5. If we created more proposals than chains, keep a log‑prob‑balanced
#      subset so that each chain explores a different mode.

dispersed_MAP_initialization <- function(
  n_inits_0,  # number of random initial points for the MAP search
  n_inits,    # number of dispersed initial points to return
  n_chains,   # chains that Stan will actually run
  stan_model, # compiled stan model (rstan::stan_model)
  stan_data,  # list of data for the model
  dirichlet_params = c(1, 1, 1, 1), # prior hyper‑parameters for theta
  alpha_params     = c(0.5, 0.05),  # mean & sd of the prior for alpha_free
  gamma_params     = c(0.1, 0.01),  # mean & sd of the prior for gamma_free
  phi_inv_lambda   = 1,             # rate of Exp(λ) prior for phi_inv_free
  beta0_params     = c(-4, 4),      # uniform bounds for beta_free
  initialize_theta   = TRUE,  # should we sample theta in step 1?
  initialize_epidem  = TRUE,  # sample epidemic‑specific params in step 1?
  jitter_init_props  = TRUE,  # jitter theta around the MAP?
  jitter_init_epidem = TRUE,  # jitter the epidemic params too?
  jitter_percentage  = 0.10,  # relative sd for Gaussian jitter
  dtv_budget         = 1,     # radius of the DTV ball used to jitter theta
  timeout           = 10 * 60, # max seconds per MAP search,
  compute_optimization = TRUE, # whether to compute the MAP
  u_map = NULL, # optional MAP vector to use instead of computing it
  parallel_threads  = 30     # max CPU cores we are willing to use
) {
  # ---- helper functions -----------------------------------------------------
  # The file below must define:
  #   * rdirichlet0()                     safe Dirichlet with zeros
  #   * run_with_timeout()                MAP optimisation with a wall‑time cap
  #   * hierarchical_DTV_ball_sampler()   exact θ‑jitter within a DTV ball
  #   * jitter_init(), jitter_init_2()    Gaussian jitter utilities
  #   * build_jittered_init()             helper to build the init list
  #   * safe_value()                       helper to extract a numeric value
  source("auxiliary functions/auxiliary_functions_initialization.R")

  available_cores <- parallel::detectCores()  # hardware concurrency

  # ---- 1. PRIOR‑BASED SEEDS -------------------------------------------------
  # Draw n_inits_0 seeds independently from the prior.  These seeds will feed
  # the short LBFGS searches that approximate the MAP.

  # θ  (mixture weights)
  theta_inits <- rdirichlet0(n_inits_0, dirichlet_params) # nolint: object_usage_linter.

  # α  (log‑latent duration E→I)
  alpha_mean <- alpha_params[1]
  alpha_sd   <- alpha_params[2]
  alpha_fixed_flag <- {
    has_fix_alpha <- "fix_alpha" %in% names(stan_data)
    has_M         <- "M" %in% names(stan_data)

    fix_alpha_val <- if (has_fix_alpha) stan_data$fix_alpha else TRUE
    M_val         <- if (has_M) stan_data$M else 0

    fix_alpha_val || M_val == 0
  }

  alpha_inits <- if (!alpha_fixed_flag) {
    rnorm(n_inits_0, alpha_mean, alpha_sd)
  } else {
    rep(0, n_inits_0)  # placeholder; Stan ignores when α is fixed
  }

  # γ  (log‑latent duration I→R)
  gamma_mean <- gamma_params[1]
  gamma_sd   <- gamma_params[2]
  gamma_fixed_flag  <- stan_data$fix_gamma
  gamma_fixed_value <- stan_data$gamma_fixed
  gamma_inits <- if (!gamma_fixed_flag) {
    rnorm(n_inits_0, gamma_mean, gamma_sd)
  } else {
    rep(gamma_fixed_value, n_inits_0)
  }

  # φ⁻¹ ~ Exp(λ)
  phi_inv_inits <- rexp(n_inits_0, phi_inv_lambda)

  # β₀  (baseline transmission)
  beta_inits <- runif(n_inits_0, beta0_params[1], beta0_params[2])

  # ---- 2. BUILD THE init LISTS ---------------------------------------------
  # rstan expects a *named list* per chain.  We assemble those lists now.
  mk_init <- function(i, with_theta, with_epi) {
    li <- list(beta_free = rep(beta_inits[i], stan_data$m - stan_data$K0))
    if (with_theta) li$theta_free <- as.numeric(theta_inits[i, ])
    if (with_epi) {
      li$alpha_free   <- alpha_inits[i]
      li$gamma_free   <- gamma_inits[i]
      li$phi_inv_free <- phi_inv_inits[i]
    }
    li
  }
  inits <- lapply(seq_len(n_inits_0), mk_init,
                  with_theta = initialize_theta,
                  with_epi   = initialize_epidem)

  # ---- 3. MAP SEARCH --------------------------------------------------------
  if(compute_optimization) {
    # Run a short optimisation from each seed.  The worker times out after 10 min
    # so that pathological seeds do not block the pool.
    optim_results <- parallel::mclapply(
      inits,
      function(init_i) {
        run_with_timeout(stan_model, stan_data, init_i, # nolint: object_usage_linter.
                        timeout = timeout, returning = "full")
      },
      mc.cores   = min(available_cores - 1, parallel_threads)
    )

    # Extract the attained log‑posterior.  Failed runs → −Inf.
    optim_logp <- vapply(optim_results,
                        safe_value, # nolint
                        numeric(1))
    MAP_index <- which.max(optim_logp)
    u_map     <- optim_results[[MAP_index]]$par  # unconstrained MAP vector
  } else if (is.null(u_map)) {
    stop("A MAP estimate was not providede or computed.")
  }

  # ---- 4. JITTER AROUND THE MAP --------------------------------------------
  # From the MAP we now craft n_inits proposals.  We have four orthogonal
  # switches controlling what gets jittered:
  #   * jitter_init_props   ⇢ θ_free via hierarchical_DTV_ball_sampler()
  #   * jitter_init_epidem  ⇢ α, γ, φ⁻¹ via Gaussian jitter
  #   * jitter_percentage   ⇢ relative width of the Gaussian jitter
  #   * dtv_budget          ⇢ radius of the DTV ball for θ

  if (jitter_init_props) {
    # Step 4A – jitter θ.
    jittered_props <- hierarchical_DTV_ball_sampler( # nolint: object_usage_linter.
      u_map$theta_free,
      dtv        = dtv_budget,
      N_Samples  = n_inits
    )
  }

  inits <- lapply(seq_len(n_inits), build_jittered_init, # nolint: object_usage_linter.
                                u_map,
                                jittered_props,
                                alpha_fixed_flag,
                                jitter_percentage,
                                gamma_fixed_flag,
                                initialize_epidem,
                                alpha_inits,
                                gamma_inits,
                                jitter_init_epidem,
                                jitter_init_props
  )

  # ---- 5. SCORE EACH PROPOSAL ----------------------------------------------
  # rstan::log_prob() needs an "empty" fit object so that it knows the model.
  fit0 <- rstan::sampling(stan_model, stan_data,
                          chains = 0, iter = 1, warmup = 0, refresh = 0)

  logps <- parallel::mclapply(
    inits,
    function(init) {
      upars <- rstan::unconstrain_pars(fit0, init)
      rstan::log_prob(fit0, upars)
    },
    mc.cores = min(available_cores - 1, n_chains)
  )

  # ---- 6. SUBSAMPLE FOR BALANCE --------------------------------------------
  # If we created more proposals than chains we cluster the log‑posterior
  # values so that each chain starts in a different high‑density basin.
  if (length(inits) > n_chains) {
    cl <- kmeans(unlist(logps), centers = n_chains)
    keep <- vapply(seq_len(n_chains),
                   function(k) sample(which(cl$cluster == k), 1),
                   integer(1))
    inits <- inits[keep]
  }

  return(list("initial_values" = inits, "logps" = unlist(logps), "MAP_params" = u_map))
}
