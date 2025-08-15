computing_and_plotting_transmission <- function(model, credible.quantiles = c(0.05, 0.95), time_basis,
                                                gamma_fixed = 0.1,
                                                params = c("p_R0", "phi_inv", "sig", "alpha",
                                                           "S0", "E0", "I0", "R0"),
                                                model_name = "SEI3R"){
  traces_model <- mcmc(get_traces(model, all = TRUE))
  traces_model <- traces_model[, params]
  # pdf(file = "synthetic_SI3R_MCMC_trace.pdf", paper = "a4", width = 0, height = 0)

  out_file <- paste0("diffusion_", model_name, "_traceplot")
  pdf(paste0(out_file, ".pdf"), width = 10, height = 14)
  par(mfrow = c(8, 2), mai = c(0.4, 0.4, 0.3, 0.3)) 
  print(plot(traces_model, auto.layout = FALSE))
  dev.off()

  out_file <- paste0("diffusion_", model_name, "_ESS")
  write.csv(effectiveSize(traces_model), file = paste0(out_file, ".csv"))
  
  posterior_model <- bi_read(model)
  
  observation_model <- sample_obs(model, start_time=Start_time, end_time=End_time, output_every=1)
  observation_test <- bi_read(observation_model)
  
  q1 <- credible.quantiles[1]
  q2 <- credible.quantiles[2]
  summary_observation_test <- summary(observation_model, type="obs", quantiles = c(q1, q2))
  # osBC_7 <- summary(obs_biBC_7, type="obs")
  summary_observation_test$time <- syn_data_difusion$time
  
  out_file <- paste0("diffusion_", model_name, "_Posterior_Predictive")
  pdf(paste0(out_file, ".pdf"), width = 14, height = 10)
  print(ggplot(summary_observation_test[1:End_time, ], aes(x = time)) +
    geom_line(aes(y = Median)) +
    geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.5) +
    geom_point(aes(y = value), syn_data_difusion[1:End_time, ], color = "darkred", size = 1) +
    ylab("Daily Incidence") + theme(axis.title.x = element_blank()))
  dev.off()
  
  beta_test <- array(dim = c(N_samples, End_time))
  R0_test <- beta_test
  
  for (np in 1:N_samples){
    beta_test[np, ] <- posterior_model$p_R0$value[np] * exp(posterior_model$log_beta$value[((np - 1)*End_time + 1):(np*End_time)])
    R0_test[np, ] <- beta_test[np, ] / gamma_fixed # posterior_model$gamma$value[np]
  }
  
  r0_test <- t(apply(data.frame(R0_test), 2, function(x) quantile(x, probs = c(0.05, 0.5, 0.95))))
  r0_test <- as.data.frame(r0_test)
  r0_test$Mean <- colMeans(data.frame(R0_test))
  r0_test$time <- 1:End_time
  r0_test$model <- exp(time_basis%*%c(-1.869874, -1.301393, -0.2422027, -1.51097, -3.304498,
                                           -3.09169, -1.568307, -1.570512, -3.447937, -4.521353,
                                           -3.334759, -2.809118)) / 0.1
  
  out_file <- paste0("diffusion_", model_name, "_Reproduction_Number")
  pdf(paste0(out_file, ".pdf"), width = 14, height = 10)
  print(ggplot(r0_test, aes(x=time)) +
    geom_line(aes(y=Mean)) +
    # geom_line(aes(y=Median)) +
    geom_ribbon(aes(ymin=`5%`, ymax=`95%`), alpha=0.5) +
    geom_line(aes(y=model), col = "green") +
    geom_hline(yintercept = 1, linetype="dashed", color = "red", size=0.5) +
    ylab("R0") + theme(axis.title.x = element_blank()))
  dev.off()
  
  beta0_test <- t(apply(data.frame(beta_test), 2, function(x) quantile(x, probs = c(q1, 0.5, q2))))
  beta0_test <- as.data.frame(beta0_test)
  beta0_test$Mean <- colMeans(beta_test)
  beta0_test$time <- 1:End_time
  beta0_test$model <- exp(time_basis%*%c(-1.869874, -1.301393, -0.2422027, -1.51097, -3.304498,
                                              -3.09169, -1.568307, -1.570512, -3.447937, -4.521353,
                                              -3.334759, -2.809118))
  
  out_file <- paste0("diffusion_", model_name, "_Beta")
  pdf(paste0(out_file, ".pdf"), width = 14, height = 10)
  print(ggplot(beta0_test, aes(x = time)) +
    geom_line(aes(y = Mean)) +
    # geom_line(aes(y=Median)) +
    geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.5) +
    geom_line(aes(y=model), col = "green") +
    ylab("beta") + theme(axis.title.x = element_blank()))
  dev.off()

  return(list(incidence = summary_observation_test, transmission = beta0_test, R0 = R0_test, trace = traces_model))
}