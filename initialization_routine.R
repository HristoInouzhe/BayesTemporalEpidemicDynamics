# Load required packages
library(rstan)
library(parallel)
library(splines)

# Set options for parallel processing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

path_to_package <- "~/Documents/" # Adjust this path to your package directory
setwd(paste0(path_to_package, "BayesTemporalEpidemicDynamics/"))

# Source the auxiliary spline functions
# (assumed to be in auxiliary_spline_functions.R)
# Folder containing your R scripts
func_dir <- "auxiliary functions"

# Find all R files in the folder
files <- list.files(func_dir, pattern = "\\.R$", full.names = TRUE)

# Source each file
invisible(lapply(files, source))

# -----------------------------
# Prepare synthetic data
# -----------------------------
seed_value <- 1234    # set a seed for reproducibility
set.seed(seed_value)  # for reproducibility

# Define the test name
test_name <- "testSEI3R_02"

data_source <- "synthetic" # "BC" for Basque Country data, "synthetic" for synthetic data

if (data_source == "BC") {
  # Load Basque Country data
  # Uncomment the following lines to load real data
  data_BC <- read.table("/data/Basque_Country_covid19_SIR_data.txt", header = FALSE)
  l_data_BC <- length(data_BC$V1)
  y <- data_BC$V1[-l_data_BC] # real observed counts
  N <- data_BC$V1[l_data_BC] # population size
} else if (data_source == "synthetic") {
  data_synthetic <- read.table("/data/synthetic_covid19_data.txt", header = FALSE)
  l_data_synthetic <- length(data_BC$V1)
  y <- data_synthetic$V1[-l_data_synthetic] # real observed counts
  N <- data_synthetic$V1[l_data_synthetic] # population size
} else {
  stop("Invalid data source specified. Use 'BC' or 'synthetic'.")
}

T_obs <- length(y)          # number of observations
ts_obs <- 1:T_obs     # observation times: 1, 2, ..., T_obs

# Model dimensions synthetic
M <- 1      # number of exposed compartments
K <- 3      # number of infected compartments
m <- 12     # number of spline basis functions
                     
# # Model dimensions BC
# M <- 1      # number of exposed compartments
# K <- 3      # number of infected compartments
# m <- 23     # number of spline basis functions

# -----------------------------
# Prepare spline objects for Stan
# -----------------------------
# Define a grid for spline evaluations;
# here we use T_obs points over the observation period
T_bs <- T_obs
ts_bs <- seq(min(ts_obs), max(ts_obs), length.out = T_bs)
bdegree <- 3       # degree for B-spline basis functions
penalty_order <- 2 # order for the random walk/difference penalty


do_under_report <- FALSE

if (do_under_report){
  under_report <- array(dim = T_obs)
  eta_0 <- 0.15;
  eta_1 <- 0.54;
  t_0 <- 92;
  t_1 <- 281;
  for (i in 1:t_0){
    under_report[i] <- eta_0
  }
  for (i in (t_0 + 1):t_1){
    under_report[i] <- eta_0 + (eta_1 - eta_0) * (i - t_0) / (t_1 - t_0) 
  }
  for (i in (t_1 + 1):T_obs){
    under_report[i] <- eta_1
  }
} else {
  under_report <- rep(1, T_obs)
}

gamma_fixed_flag <- TRUE
gamma_fixed_value <- 1 / 10 # 1 / (mean infectious period)

n_chains <- 10 # number of chains
n_inits_0 <- 10 * n_chains # number of initializations
n_inits <- 10 * n_chains

a_1 <- 999993.424608; a_2 <- 4.575392; a_3 <- 1.0; a_4 <- 1.0 # For synthetic data a_0 = 1e+6, E[E_0] = 10
alpha_mean <- 0.5; alpha_sd <- 0.05
phi_inv_lambda <- 20

initialize_theta <- TRUE # whether to initialize theta_free
initialize_epidem <- TRUE # whether to initialize epidemic parameters

MAP_initialization <- TRUE # use MAP initialization

jitter_percentage <- 0.25 # percentage of the initial theta to jitter
dtv_budget <- 1e-4 # radius of the DTV ball used to jitter theta
jitter_init_props <- TRUE # whether to jitter the initial proportions
jitter_init_epidem <- TRUE # whether to jitter the initial epidemic parameters

K0 <- 0 # number of first beta's to fix (0 means no fixing)

timeout <- 10 * 60 # max seconds per MAP search
parallel_threads <- 20 # number of threads for parallel processing

save_to_haics <- FALSE # whether to save the initialization results for HAICS

# Use the auxiliary function to compute
# the B-spline basis matrix and penalty matrix
spline_objs <- compute_spline_objects(ts_bs, bdegree, m, penalty_order)

# -----------------------------
# Prepare the data list for Stan
# -----------------------------

stan_data <- list(
  # ─── observation data ───
  T_obs         = T_obs,
  t0            = 0,
  ts_obs        = as.numeric(ts_obs),      # numeric vector length T_obs
  y             = as.integer(y),           # integer vector length T_obs
  eta           = under_report,            # under report level

  # ─── model dims & splines ───
  M             = M,            # integer
  K             = K,            # integer
  m             = m,            # integer (# basis fns)
  T_bs          = T_bs,         # integer
  ts_bs         = as.numeric(spline_objs$ts_bs),  # numeric vector length T_bs
  B             = spline_objs$B_matrix,           # matrix T_bs × m
  K_bs          = spline_objs$K,                  # matrix m × m

  # ─── priors & constants ───
  N             = N,            # pop size
  alpha_dir     = c(a_1, a_2, a_3, a_4),    # numeric length 4
  prior_mean_alpha   = alpha_mean,  # numeric
  prior_sd_alpha     = alpha_sd, # numeric > 0
  prior_mean_gamma   = 0.1,  # numeric
  prior_sd_gamma     = 0.01, # numeric > 0
  prior_rate_phi_inv = phi_inv_lambda, # numeric > 0

  # ─── fix switches for (α, γ, φ⁻¹, β) ───
  fix_alpha     = 0,            # 0 = estimate α, 1 = fix
  alpha_fixed   = 0.5,        # whatever you like
  fix_gamma     = 1,            # 0 = estimate γ, 1 = fix
  gamma_fixed   = gamma_fixed_value,
  fix_phi_inv   = 0,
  phi_inv_fixed = 0.1,
  K0            = K0,            # how many of first β’s to fix

  # ─── fix switches for initial states ───
  fix_S0        = 0,            # whether to fix S0
  S0_fixed      = (N - 10) / N, # value of S0 if fixed
  fix_E0        = 0,            # whether to fix E0
  E0_fixed      = 10 / N,       # value of E0 if fixed
  fix_I0        = 0,            # whether to fix I0
  I0_fixed      = 0             # valeue of I0 if fixed
)

# Build a dummy m‐vector, then overwrite the first K0 elements
stan_data$beta_fixed <- rep(0, stan_data$m)
if (K0 > 0) {
  stan_data$beta_fixed[1:K0] <- c(-1.8699, -1.3014)[1:K0]
}

# -----------------------------
# Compile the Stan model
# -----------------------------
stan_model <- stan_model(file = "Stan/model/SEMIKR_Dirichlet_final.stan") # For the model with E compartment
# stan_model <- stan_model(file = "Stan/model/SIKR_Dirichlet_final.stan") # For the model without E compartment


if (MAP_initialization{
  inits <- dispersed_MAP_initialization(
  n_inits_0          = n_inits_0,
  n_inits            = n_inits,
  n_chains           = n_chains,
  stan_model         = stan_model,
  stan_data          = stan_data,
  dirichlet_params   = c(a_1, a_2, a_3, a_4), # prior hyper-parameters for theta, must be length 4 for model with E compartment, and 3 for model without E
  alpha_params       = c(alpha_mean, alpha_sd),
  gamma_params       = c(0.1, 0.01),
  phi_inv_lambda     = phi_inv_lambda,
  beta0_params       = c(-4, 4),
  jitter_percentage  = jitter_percentage,
  dtv_budget         = dtv_budget,
  timeout           = timeout, # max seconds per MAP search
  compute_optimization = TRUE, # whether to compute the MAP
  parallel_threads   = parallel_threads #parallel_threads
)
})

save(inits, file = paste0("initializations/", test_name, "_inits.RData")) # Save the initialization results
save(stan_data, file = paste0("Stan/stan_data_", test_name, ".RData")) # Save the Stan data

if(save_to_haics) {
  # Save the initialization results for HAICS
  input_file_path <- paste0("/haics/simulation/output/inputfile_", test_name, ".txt") # must be defined before calling this function and in the correct HAICS format
  saving_initialization_to_haics(inits, N, test_name, input_file_path)
}
