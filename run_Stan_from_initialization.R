library(rstan)

stan_model <- stan_model(file = "Stan/model/SEMIKR_Dirichlet_final.stan") # For the model with E compartment
# stan_model <- stan_model(file = "Stan/model/SIKR_Dirichlet_final.stan") # For the model without E compartment

experiment_id <- "testSEI3R_02" # or "synthetic"
load(paste0("initializations/", experiment_id, "_inits.RData")) # Load the initialization results
load(paste0("Stan/stan_data_", experiment_id, ".RData")) # Load the Stan data

n_chains <- 4 # number of chains
iter <- 2000 # total iterations per chain: 1000 warmup + 1000 sampling
warmup <- 1000 # number of warmup iterations per chain
seed_value <- 1234 # for reproducibility

fit <- sampling(
  stan_model,       # compiled Stan model
  data = stan_data, # data list prepared above
  chains = n_chains,      # number of chains
  init = inits$initial_values,     # initial values for each chain
  iter = iter,      # total iterations per chain: 1000 warmup + 1000 sampling
  warmup = warmup,  # number of warmup iterations per chain
  seed = seed_value,      # for reproducibility
  refresh = 50,    # progress output every 100 iterations
  control = list(
    adapt_delta = 0.75, #0.9
    # target acceptance rate, default is 0.8,
    # 0.65 for exploring faster, 0.95 for more accurate
    adapt_init_buffer = 150, #500
    # increase initial fast‐adapt interval, default is 75
    adapt_term_buffer = 100,
    # decrease final fast‐adapt interval, default is 50
    max_treedepth = 15, # maximum tree depth, default is 10
    stepsize = 1,   # initial step size, default is 1, # 2
    stepsize_jitter = 0.35
    # jitter for step size, default is 0.05
    # [ε×(1−stepsize_jitter),ε×(1+stepsize_jitter)]
  ),
  verbose = TRUE    # show progress
)

# Save the fitted model
saveRDS(fit, file = paste0("Stan/stan_fit_", experiment_id, ".RDS"))

plot_Stan_fit(paste0("Stan/stan_fit_", experiment_id, ".RDS"), stan_data$y)