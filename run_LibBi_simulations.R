library(rbi)
library(coda)
library(ggplot2)
library(MASS)
library(Matrix)
library(splines)
library(readr)

source("auxiliary functions/posteriorSamplingRbi.R")
source("auxiliary functions/computing_and_plotting_transmission.R")
source("auxiliary functions/auxiliary_spline_functions.R")

N_samples = 1000
Start_time = 1
End_time = 100
N_particles = 1000

data_synthetic <- read.table("/data/synthetic_covid19_data.txt", header = FALSE)
l_data_synthetic <- length(data_BC$V1)
y <- data_synthetic$V1[-l_data_synthetic] # real observed counts

syn_data_difusion <- data.frame(Incidence = y,
                                 Time = Start_time:End_time)

# the B-spline basis matrix and penalty matrix
spline_objs <- compute_spline_objects(1:End_time, 3, 12, 2)
time_basis <- spline_objs$B_matrix

model_test_SIR <- posteriorSamplingRbi("SIR_beta_smooth_stan_3.bi",
                                       list(Incidence = syn_data_difusion), 
                                       "SIR_synthetic",
                                       data.frame(N_samples = N_samples,
                                                  N_particles = N_particles,
                                                  End_time = End_time))


posterior_prediction_SIR <- computing_and_plotting_transmission(model_test_SIR, time_basis = time_basis,
                                                                   gamma_fixed = 0.1,
                                                                   params = c("p_R0", "phi_inv", "sig", "alpha",
                                                                              "S0", "I0", "R0"),
                                                                   model_name = "SIR")


model_test_SI3R <- posteriorSamplingRbi("SI3R_beta_smooth_stan_3.bi",
                                        list(Incidence = syn_data_difusion), 
                                        "SI3R_synthetic",
                                        data.frame(N_samples = N_samples,
                                                   N_particles = N_particles,
                                                   End_time = End_time))

posterior_prediction_SI3R <- computing_and_plotting_transmission(model_test_SI3R_v2, time_basis = time_basis,
                                                                    gamma_fixed = 0.1,
                                                                    params = c("p_R0", "phi_inv", "sig", "alpha",
                                                                               "S0", "I0", "R0"),
                                                                    model_name = "SI3R")

model_test_SEIR <- posteriorSamplingRbi("SEIR_beta_smooth_stan_3.bi",
                                        list(Incidence = syn_data_difusion), 
                                        "SEIR_synthetic",
                                        data.frame(N_samples = N_samples,
                                                   N_particles = N_particles,
                                                   End_time = End_time))

posterior_prediction_SEIR <- computing_and_plotting_transmission(model_test_SEIR_v2, time_basis = time_basis,
                                                                 gamma_fixed = 0.1,
                                                                 params = c("p_R0", "phi_inv", "sig", "alpha",
                                                                            "S0", "E0", "I0", "R0"),
                                                                 model_name = "SEIR")

model_test_SEI3R <- posteriorSamplingRbi("SEI3R_beta_smooth_stan_3.bi",
                                         list(Incidence = syn_data_difusion), 
                                         "SEI3R_synthetic",
                                         data.frame(N_samples = N_samples,
                                                    N_particles = N_particles,
                                                    End_time = End_time))

posterior_prediction_SEI3R <- computing_and_plotting_transmission(model_test_SEI3R, time_basis = time_basis,
                                                                     gamma_fixed = 0.1,
                                                                     params = c("p_R0", "phi_inv", "sig", "alpha",
                                                                                         "S0", "E0", "I0", "R0"),
                                                                     model_name = "SEI3R")