posteriorSamplingRbi <- function(model_bi, data, model_name, sampling_parameters){
  N_samples <- sampling_parameters$N_samples
  N_particles <- sampling_parameters$N_particles
  End_time <- sampling_parameters$End_time
  
  SIRmodel_daily_test <- bi_model(model_bi)
  model_test <- libbi(SIRmodel_daily_test)
  
  t0 <- Sys.time()
  model_test <- sample(model_test, target="posterior", nsamples = N_samples,
                           nparticles = N_particles, end_time = End_time,
                           obs = data,
                           nthreads = 8, verbose = TRUE, sampler = "sir") # 
  t1 <- Sys.time()
  t_test <- t1 - t0
  print(t_test)
  save_libbi(model_test, name = model_name)
  return(model_test)
}
  