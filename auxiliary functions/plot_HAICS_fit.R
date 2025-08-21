plot_HAICS_fit <- function(fit_file, y) {
  require(outbreaks)
  require(ggplot2)
  require(splines)
  require(coda)
  require(invgamma)
  require(RColorBrewer)
  require(parallel)
  require(data.table)
  require(MASS)


  source("analysis_haics_results.R")
  source("auxiliary_ploting_hmc.R")
  source("auxiliary_spline_functions.R")
  load("medidas_cierre_apertura.RData")
  load("syn_data_difusion.RData")
  load("synthetic_data_erlang_exposed.RData")
  load("Incidence_Data.RData")


  experiment_name <- "testSens_07"
  inptu_file_data <- read_haics_config(paste0("/nvme0n1-disk/haics/simulation/output/inputfile_", experiment_name,".txt"))
  pipeline_params <- set_pipeline_params(inptu_file_data, chains_default = 10,
                                        file_directory_default = "/nvme0n1-disk/haics/simulation/output/",
                                        file_name_default = experiment_name,
                                        npop_default = 2185605)
  t_fin <- pipeline_params$t_fin
  time_basis_0 <- pipeline_params$time_basis_0
  n_i_knots <-pipeline_params$n_i_knots
  time_basis <- pipeline_params$time_basis
  Gamma.fixed <- pipeline_params$Gamma.fixed
  chains <- pipeline_params$chains
  file_directory <- pipeline_params$file_directory
  model_type <- pipeline_params$model_type
  dirich_alpha <- pipeline_params$dirich_alpha
  burnin <- pipeline_params$burnin
  gamma.restrictions <- pipeline_params$gamma.restrictions
  n_var <- pipeline_params$n_var
  neq <- pipeline_params$neq
  plotTitle <- pipeline_params$plotTitle
  syn_data <- pipeline_params$syn_data 
  ndata <- pipeline_params$ndata 
  npop <- pipeline_params$npop
  do_plot <- FALSE
  measures <- medidas_cierre_apertura
  # Determine the number of cores to use (leaving one free for the system)
  no_cores <- min(max(detectCores() - 1, 1), chains)
  # no_cores <- 3

  ################################################################################
  ################################################################################

  if (gamma.restrictions[1] == 1) {
    # Gamma is fixed with value gamma.restrictions[2]
    Gamma.fixed <- gamma.restrictions[2]
    # Determine the number of cores to use (leaving one free for the system)
    # no_cores <- min(max(detectCores() - 1, 1), chains)
    
    # Parallel processing of chains using mclapply
    results <- mclapply(1:chains, function(i) {
      traj_data <- tryCatch(
        {
          fread(paste0(file_directory, file_name, "_", i, "/trajectories.txt"))
        },
        error = function(e) {
          message("An error occurred while reading the file: ", e$message)
          NULL  # Return NULL or any default value you prefer
        }
      )
      ode_data <- tryCatch(
        {
          fread(paste0(file_directory, file_name,"_", i,"/ode_sol", ".txt"))
        },
        error = function(e) {
          message("An error occurred while reading the file: ", e$message)
          NULL  # Return NULL or any default value you prefer
        }
      )
        
      if (is.null(traj_data) || is.null(ode_data)){
        list("Error")
      } else{
        n_traj = dim(traj_data)[1]
        trajectories_original <- if(model_type == "SIR"){
          data.frame(
            Gamma = traj_data[, 1],
            S0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 2])),
            I0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 3])),
            Phi_Inv = exp(traj_data[, 4]) + 1e-10,
            Tau = exp(traj_data[, 5]),
            traj_data[, 6:n_var]
          )
        } else {
          data.frame(
            Alpha = exp(traj_data[, 1]) + 1e-10,
            Gamma = traj_data[, 2],
            S0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 3])),
            E0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 4])),
            I0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 5])),
            Phi_Inv = exp(traj_data[, 6]) + 1e-10,
            Tau = exp(traj_data[, 7]),
            traj_data[, 8:n_var]
          )
        }
        
        list(trajectories_original = trajectories_original, ode_sol = ode_data)
      }
    }, mc.cores = no_cores)
  } else {
    if (gamma.restrictions[3] == 1){
      # gamma is bounded with bounds given in gamma_bounds
      # Defining the bounds for the gamma parameter 
      gammaL <- gamma_bounds[1]
      gammaU <- gamma_bounds[2]
      
      # Determine the number of cores to use (leaving one free for the system)
      # no_cores <- min(max(detectCores() - 1, 1), chains)
      
      # Parallel processing of chains using mclapply
      results <- mclapply(1:chains, function(i) {
        traj_data <- tryCatch(
          {
            fread(paste0(file_directory, file_name, "_", i, "/trajectories.txt"))
          },
          error = function(e) {
            message("An error occurred while reading the file: ", e$message)
            NULL  # Return NULL or any default value you prefer
          }
        )
        ode_data <- tryCatch(
          {
            fread(paste0(file_directory, file_name,"_", i,"/ode_sol", ".txt"))
          },
          error = function(e) {
            message("An error occurred while reading the file: ", e$message)
            NULL  # Return NULL or any default value you prefer
          }
        )
        if (is.null(traj_data) || is.null(ode_data)){
          list("Error")
        } else{
          n_traj = dim(traj_data)[1]
          trajectories_original <- if(model_type == "SIR"){
            data.frame(
              S0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 2])),
              I0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 3])),
              Phi_Inv = exp(traj_data[, 4]) + 1e-10,
              Tau = exp(traj_data[, 5]),
              traj_data[, 6:n_var]
            )
          } else {
            data.frame(
              Alpha = exp(traj_data[, 1]) + 1e-10,
              Gamma = gammaL + (gammaU - gammaL) / (1 + exp(-traj_data[, 2])),
              S0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 3])),
              E0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 4])),
              I0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 5])),
              Phi_Inv = exp(traj_data[, 6]) + 1e-10,
              Tau = exp(traj_data[, 7]),
              traj_data[, 8:n_var]
            )
          }
          list(trajectories_original = trajectories_original, ode_sol = ode_data)
        }
      }, mc.cores = no_cores)
    } else{
      # gamma is not fixed or bounded
      
      # Determine the number of cores to use (leaving one free for the system)
      # no_cores <- min(max(detectCores() - 1, 1), chains)
      
      # Parallel processing of chains using mclapply
      results <- mclapply(1:chains, function(i) {
        traj_data <- tryCatch(
          {
            fread(paste0(file_directory, "/", file_name, "_", i, "/trajectories.txt"))
          },
          error = function(e) {
            message("An error occurred while reading the file: ", e$message)
            NULL  # Return NULL or any default value you prefer
          }
        )
        ode_data <- tryCatch(
          {
            fread(paste0(file_directory, file_name, "_", i, "/ode_sol", ".txt"))
          },
          error = function(e) {
            message("An error occurred while reading the file: ", e$message)
            NULL  # Return NULL or any default value you prefer
          }
        )
        
        if (is.null(traj_data) || is.null(ode_data)){
          list("Error")
        } else{
          n_traj = dim(traj_data)[1]
          trajectories_original <- if(model_type == "SIR"){
            data.frame(
              Gamma = exp(traj_data[, 1]) + 1e-10,
              S0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 2])),
              I0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 3])),
              Phi_Inv = exp(traj_data[, 4]) + 1e-10,
              Tau = exp(traj_data[, 5]),
              traj_data[, 6:n_var]
            )
          } else {
            data.frame(
              Alpha = exp(traj_data[, 1]) + 1e-10,
              Gamma = exp(traj_data[, 2]) + 1e-10,
              S0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 3])),
              E0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 4])),
              I0 = 0 + (npop - 0) / (1 + exp(-traj_data[, 5])),
              Phi_Inv = exp(traj_data[, 6]) + 1e-10,
              Tau = exp(traj_data[, 7]),
              traj_data[, 8:n_var]
            )
          }
          list(trajectories_original = trajectories_original, ode_sol = ode_data)
        }
      }, mc.cores = no_cores)
    }
  }

  # Unpack results
  trajectories_original <- lapply(results, `[[`, "trajectories_original")
  ode_sol <- lapply(results, `[[`, "ode_sol")

  bad_chains <- c() 
  if(length(bad_chains) > 0){
    valid_chains <- (1:chains)[-bad_chains]
  } else {
    valid_chains <- (1:chains)
  }


# Initialize an empty list to store initial splines and parameters
  initial.splines <- list() 
  initial.params <- list() 
  ii <- 0
  for (i in valid_chains){
    ii <- ii + 1
    if(model_type == "SEIR"){
      # Read and store the data from a file containing initial splines
      initial.splines[[ii]] <- t(read.table(paste0(file_directory, file_name, "_", i, "/splinebasis_", file_name, "_", i, ".txt"), quote="\"", comment.char=""))

      
      # Read and store the data from a file containing initial parameters
      initial.params.temp <- read.table(paste0(file_directory, file_name, "_", i, "/initialpoint_", file_name, "_", i, ".txt"), quote="\"", comment.char="")
      names(initial.params.temp) <- c("alpha", "gamma", "S0", "E0", "I0", "phi_inv", "tau")
      initial.params[[ii]] <- initial.params.temp
    } else {
      # Read and store the data from a file containing initial splines
      initial.splines[[ii]] <- t(read.table(paste0(file_directory, file_name, "_", i, "/splinebasis_", file_name, "_", i, ".txt"), quote="\"", comment.char=""))
      
      # Read and store the data from a file containing initial parameters
      initial.params.temp <- read.table(paste0(file_directory, file_name, "_", i, "/initialpoint_", file_name, "_", i, ".txt"), quote="\"", comment.char="")
      names(initial.params.temp) <- c("gamma", "S0", "I0", "phi_inv", "tau")
      initial.params[[ii]] <- initial.params.temp
    }
  }
  remove(initial.params.temp)
  remove(results)
  gc()

  # Initialize an array to store the Effective Sample Size (ESS) for each chain and each variable.
  # 'chains' is the number of MCMC chains, and 'n_var' is the number of variables.
  if (gamma.restrictions[1] == 1){
    ESS_chains <- array(dim = c(chains -length(bad_chains), n_var - 1)) # gamma is fixed
  } else{
    ESS_chains <- array(dim = c(chains -length(bad_chains), n_var)) # gamma not fixed
  }
  hmc_list <- list()

  # Loop over each chain.
  ii <- 0
  for (i in valid_chains){
    ii = ii + 1
    # Convert the trajectory of the i-th chain to an MCMC object, excluding the burn-in period.
    # 'burn_in' is the number of initial steps of the chain to be discarded.
    # 'trajectories_original[[i]]' contains the trajectory of the i-th chain.
    if (gamma.restrictions[1] == 1){
      mcmc_obj <- as.mcmc(trajectories_original[[i]][(burnin+1):dim(trajectories_original[[i]])[1], -2])
    } else{
      mcmc_obj <- as.mcmc(trajectories_original[[i]][(burnin+1):dim(trajectories_original[[i]])[1], ])
    }
    
    # Calculate the Effective Sample Size for each variable in the i-th chain.
    # The 'effectiveSize' function is part of the 'coda' package and calculates ESS.
    ESS_chains[ii, ] <- effectiveSize(mcmc_obj)
    hmc_list[[ii]] <- mcmc_obj
  }
  R_hat_values <- gelman.diag(as.mcmc.list(hmc_list), transform = TRUE, autoburnin = FALSE)
  remove(hmc_list)
  print(ESS_chains)
  print(R_hat_values)
  if (model_type == "SEIR"){
    write.csv(R_hat_values$psrf, file = paste0(file_directory, file_name, "_R_hat.csv"),
              row.names = c("alpha", "S0", "E0", "I0", "Phi_Inv", "Tau",
                            paste0(rep("beta_", length(time_basis[1, ])), 1:length(time_basis[1, ]))))
    colnames(ESS_chains) <- c("alpha", "S0", "E0", "I0", "Phi_Inv", "Tau",
                              paste0(rep("beta_", length(time_basis[1, ])), 1:length(time_basis[1, ])))
    write.csv(ESS_chains, file = paste0(file_directory, file_name, "_ESS.csv"))
  } else if (model_type == "SIR"){
    write.csv(R_hat_values$psrf, file = paste0(file_directory, file_name, "_R_hat.csv"),
              row.names = c("S0", "I0", "Phi_Inv", "Tau",
                            paste0(rep("beta_", length(time_basis[1, ])), 1:length(time_basis[1, ]))))
    colnames(ESS_chains) <- c("S0", "I0", "Phi_Inv", "Tau",
                              paste0(rep("beta_", length(time_basis[1, ])), 1:length(time_basis[1, ])))
    write.csv(ESS_chains, file = paste0(file_directory, file_name, "_ESS.csv"))
  }


  # Create empty lists for accepted proposals and ODE solutions
  accepted_proposals <- list()
  accepted_ode_sol <- list()

  # Process each chain for accepted proposals
  ii <- 0
  for (i in valid_chains) {
    ii <- ii + 1
    start_idx = ifelse(burnin == 0, 1, burnin + 1)
    
    # Extract the relevant part of trajectories_original after burn-in
    relevant_traj = trajectories_original[[i]][start_idx:nrow(trajectories_original[[i]]), ]
    
    # Get unique row indices from the relevant part of trajectories_original
    unique_indices = which(!duplicated(relevant_traj))
    
    # Adjust indices to align with original ode_sol (since unique_indices are based on post-burn-in subset)
    adjusted_indices = unique_indices + burnin - 1
    
    # Use the adjusted unique indices to extract rows from trajectories_original and ode_sol
    accepted_proposals[[ii]] <- relevant_traj[unique_indices, ]
    accepted_ode_sol[[ii]] <- ode_sol[[i]][adjusted_indices, ]
  }

  # Create an array to store the dimensions of accepted proposals for each chain
  dim_a_prop <- array(0, dim = valid_chains)

  # Create an array to store the acceptance rate for each chain
  acceptance_rate <- array(0, dim = valid_chains)

  # Loop over each chain
  ii <- 0
  for (i in valid_chains) {
    ii <- ii + 1
    # Calculate the number of rows (dimensions) of accepted proposals for the current chain
    dim_a_prop[ii] <- dim(accepted_proposals[[ii]])[1]
    
    # Calculate the acceptance rate for the current chain
    acceptance_rate[ii] <- dim_a_prop[ii] / (dim(trajectories_original[[i]])[1] - burnin)
  }

  # Check if the number of chains is 10 or less.
  if (length(valid_chains) <= 10){
    # If there are 10 or fewer chains, set a color palette using the 'brewer.pal' function.
    # The palette is set to have 10 colors from the "Paired" palette.
    palette(brewer.pal(n = 10, name = "Paired"))
  } else {
    # If there are more than 10 chains, set a color palette with a number of colors
    # equal to the number of chains. This ensures each chain gets a unique color.
    palette(brewer.pal(n = length(valid_chains), name = "Paired"))
  }

  # Initialize 'fills_plot' with the first color from the palette set above.
  fills_plot <- c(palette()[1])

  # Check if the number of 'chains' is greater than 1.
  if(chains > 1){
    # If there are multiple chains, loop through each chain starting from the second.
    for(i in 2:chains){
      # For each chain, add the corresponding color from the palette to 'fills_plot'.
      # This loop collects colors for each chain to be used in the plot.
      fills_plot <- c(fills_plot, palette()[i])
    } 
  }

  # Initialize a vector with the name for the first chain.
  colnames_fills_plot <- c("Chain_1")

  # Check if there are more than one chain.
  if(chains > 1){
    # If there are multiple chains, create a label for each chain.
    for (i in 2:chains){
      # Append the name for each subsequent chain (Chain_2, Chain_3, etc.) to the vector.
      colnames_fills_plot <- c(colnames_fills_plot, paste0("Chain_", i))
    }
  }

  # Assign these names to 'fills_plot'. 
  names(fills_plot) <- colnames_fills_plot

  print("Starting Posterior predictive check")
  t0 <- Sys.time()
  # Calculate the incidence based on the solutions from the differential equations.
  # 'accepted_ode_sol' contains the solutions from ODE (Ordinary Differential Equation) models for each chain.
  # This calculation is done in parallel using multiple cores for efficiency.
  model_incidence <- parallel::mclapply(accepted_ode_sol, function(y){
    as.data.frame(t(apply(y, 1, function(x) {
      x[seq(neq, ndata * neq, neq)] - c(0, x[seq(neq, (ndata - 1) * neq, neq)])
    })))
  }, mc.cores = no_cores)

  # Perform posterior predictive checks.
  # This process differs based on the model type (SIR or another type).
  posterior_predictive_check <- mclapply(1:length(valid_chains), function(j) {
    posterior <- array(dim = c(dim_a_prop[[j]], ndata))
    
    for (i in 1:dim_a_prop[[j]]) {
      # Determine the 'size' parameter based on the model type
      size_param = if(model_type == "SIR") trajectories_original[[valid_chains[j]]][i, 3 + 1] else trajectories_original[[valid_chains[j]]][i, 4 + 2]
      
      # Vectorized generation of random numbers using rnegbin
      posterior[i, ] <- rnegbin(as.numeric(model_incidence[[j]][i, ]), theta = 1 / size_param)
    }
    return(posterior)
  }, mc.cores = no_cores)

  # Initialize a data frame for the first chain's posterior predictive check results.
  j <- 1 # valid_chains[1]
  provisional <- cbind(
    # Calculate the 5th, 50th (median), and 95th percentiles and the mean for each time point.
    as.data.frame(t(apply(posterior_predictive_check[[j]], 2, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))),
    colMeans(posterior_predictive_check[[j]], na.rm = TRUE)
  )
  colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"))
  posterior_predictive_check_plot2 <- provisional

  for (j in 2:length(valid_chains)){ # valid_chains[-1]
    provisional <- cbind(
      as.data.frame(t(apply(posterior_predictive_check[[j]], 2, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))),
      colMeans(posterior_predictive_check[[j]], na.rm = TRUE)
    )
    colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"))
    posterior_predictive_check_plot2 <- cbind(posterior_predictive_check_plot2, provisional)
  }

  # Add observed data to the data frame.
  if(syn_data){
    # If the data is synthetic, include the synthetic data for comparison.
    posterior_predictive_check_plot2 <- cbind(
      data.frame(
        day = 1:ndata, 
        data = synthetic_data_erlang_exposed$Daily_Incidence_NB, 
        Mean = synthetic_data_erlang_exposed$Daily_Incidence
      ), 
      posterior_predictive_check_plot2
    )
  } else {
    # If the data is not synthetic, include the actual observed data for comparison.
    posterior_predictive_check_plot2 <- cbind(
      data.frame(
        day = Incidence_Data$Date, 
        data_corr = Incidence_Data$Modif_Cases[1:ndata], 
        data = Incidence_Data$Cases[1:ndata]
      ), 
      posterior_predictive_check_plot2
    )
  }

  # Initialize the total posterior predictive check data with the first chain's data.
  posterior_predictive_check_total <- posterior_predictive_check[[1]] # valid_chains[1]

  # Loop over the remaining chains to append their posterior predictive check data to the total data.
  for(iii in 2:length(valid_chains)){ # valid_chains[-1]
    posterior_predictive_check_total <- rbind(posterior_predictive_check_total, posterior_predictive_check[[iii]])
  }

  # Conditional plotting based on whether the data is synthetic.
  if(syn_data){
    # Initialize the ggplot object for synthetic data.
    graph_to_plot_3 <- ggplot(posterior_predictive_check_plot2, aes(x = day))
    
    # Loop over the chains to add graphical elements for each chain.
    for(i in 1:length(valid_chains)){ #valid_chains # 1:chains
      # Add a ribbon for the 5th to 95th percentile range and a line for the median.
      graph_to_plot_3 <- graph_to_plot_3 + 
        geom_ribbon(aes(ymin = !!sym(paste0("Chain", i, "_5")),
                        ymax = !!sym(paste0("Chain", i, "_95"))), fill = fills_plot[i], alpha = 0.4) +
        geom_line(aes(y = !!sym(paste0("Chain", i, "_Mean"))), color = fills_plot[i], linewidth = 1.15)
    }
    
    # Add more elements to the plot.
    graph_to_plot_3 <- graph_to_plot_3 + 
      geom_point(aes(y = data)) +
      geom_line(aes(y = Mean), color = "green", linetype = 2, linewidth = 1.5) +
      ylab("Daily Incidence") + theme(axis.title.x = element_blank(), legend.position = "none") #+
      #scale_fill_manual(values = fills_plot) + scale_color_manual(values = fills_plot)
  } else {
    # Initialize the ggplot object for non-synthetic data.
    graph_to_plot_3 <- ggplot(posterior_predictive_check_plot2, aes(x = day))
    
    # Similar loop for adding elements for each chain.
    for(i in 1:length(valid_chains)){
      graph_to_plot_3 <- graph_to_plot_3 + 
        geom_ribbon(aes(ymin = !!sym(paste0("Chain", i, "_5")),
                        ymax = !!sym(paste0("Chain", i, "_95"))), fill = fills_plot[i], alpha = 0.5) +
        geom_line(aes(y = !!sym(paste0("Chain", i, "_50"))), color = fills_plot[i], linewidth = 1.15)
    }
    
    # Add more elements to the plot, including points for observed data and vertical lines.
    graph_to_plot_3 <- graph_to_plot_3 + 
      geom_point(aes(y = data_corr)) +
      geom_point(aes(y = data), color = "orange") +
      geom_vline(xintercept = as.Date(measures$Date), linetype = "dotdash",
                color = measures$Colour, linewidth = 0.5) +
      ylab("Daily Incidence") + xlab("Date") + theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
                                                    axis.text = element_text(size = 20),
                                                    axis.title = element_text(size = 27), legend.position = "none") +
      labs(title = paste0("Posterior predictive check: ", plotTitle)) +
      scale_fill_manual(values = fills_plot) + coord_cartesian(ylim = c(0, 6000)) +
      scale_x_date(breaks = "month", date_labels = "%m/%y") + scale_color_manual(values = fills_plot)
  }

  # Print the plot if enabled.
  if(do_plot){
    print(graph_to_plot_3)
  }
  pdf(file = paste0(file_directory,"plots/",file_name, "_posteriorpredictive.pdf"), paper = "a4r", width = 0, height = 0)
  print(graph_to_plot_3)
  dev.off()

  t1 <- Sys.time()
  print(paste0("Finished Posterior Predictive Check: ", t1 - t0))

  remove(ode_sol)
  remove(accepted_ode_sol)
  remove(posterior_predictive_check_total)
  remove(posterior_predictive_check_plot2)
  remove(posterior_predictive_check)
  remove(model_incidence)
  remove(graph_to_plot_3)
  remove(provisional)
  remove(relevant_traj)
  gc()

  print("Starting Beta computation")
  t0 <- Sys.time()
  # Initialize an empty list called 'beta' (the time-dependent transmission rate).
  beta <- list()

  # Check the type of model being used.
  if(model_type == "SIR"){
    # If the model type is "SIR", perform a calculation for each chain.
    for (j in 1:length(valid_chains)) { # 1:chains
      # Calculate 'beta' for each chain using a matrix product ('%*%') between 'time_basis'
      # and a subset of 'accepted_proposals', then exponentiate the result.
      beta[[j]] <- exp(time_basis %*% t(accepted_proposals[[j]][, 6:n_var]))
    }
  } else {
    # If the model type is not "SIR", use a similar calculation but start from the 6th column.
    for (j in 1:length(valid_chains)) {# 1:chains
      beta[[j]] <- exp(time_basis %*% t(accepted_proposals[[j]][, 8:n_var]))
    }
  }

  # Initialize a variable for the first chain.
  j <- 1 #valid_chains[1]

  # Calculate quantiles and mean for the first chain and store it in a data frame.
  # The 'apply' function is used to calculate quantiles and mean for each row in 'beta[[j]]'.
  provisional <- cbind(as.data.frame(t(apply(beta[[j]], 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))), 
                      rowMeans(beta[[j]], na.rm = TRUE), exp(time_basis %*% initial.splines[[j]]))
  # Set column names for the 'provisional' dataframe.
  colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"), paste0("starting_", j))

  # Initialize 'beta_plot2' (which will be used for plotting the transmission rate)
  # with the values from 'provisional'.
  beta_plot2 <- provisional

  # Process data for remaining chains.
  for (j in 2:length(valid_chains)){ # valid_chains[-1]
    provisional <- cbind(as.data.frame(t(apply(beta[[j]], 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))), 
                        rowMeans(beta[[j]], na.rm = TRUE), exp(time_basis %*% initial.splines[[j]]))
    colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"), paste0("starting_", j))
    beta_plot2 <- cbind(beta_plot2, provisional)
  }

  # Conditional plotting based on whether the data is synthetic ('syn_data') or not.
  if (syn_data){
    # If the data is synthetic, bind additional columns (the true transmission rate) to 'beta_plot2'.
    beta_plot2 <- cbind(data.frame(
      day = 1:ndata, 
      data = exp(time_basis_0 %*% c(-1.869874, -1.301393, -0.2422027, -1.51097, -3.304498, -3.09169, -1.568307, -1.570512, -3.447937, -4.521353, -3.334759, -2.809118))
    ), beta_plot2)
  } else {
    # If the data is not synthetic, bind different columns (the initial proposal for the MCMC) to 'beta_plot2'.
    beta_plot2 <- cbind(data.frame(
      day = Incidence_Data$Date
    ), beta_plot2)
  }

  # Combine all chains into one data structure.
  beta_total <- beta[[1]]
  for(iii in 2:length(valid_chains)){ # 2:chains
    beta_total <- cbind(beta_total, beta[[iii]])
  }

  # Conditional plot generation based on whether the data is synthetic or not.
  if(syn_data){
    # Start building a plot with 'day' on the x-axis.
    graph_to_plot <- ggplot(beta_plot2, aes(x = day))
    
    # Loop over the number of chains to add graphical elements for each chain.
    for(i in 1:length(valid_chains)){ #valid_chains
      # Add a ribbon (shaded area) for the 5th and 95th percentiles.
      graph_to_plot <- graph_to_plot + geom_ribbon(aes(ymin = !!sym(paste0("Chain", i, "_5")),
                                                      ymax = !!sym(paste0("Chain", i, "_95"))),
                                                  fill = fills_plot[i], group =  factor(paste0("Chain_", i)),
                                                  alpha = 0.4) + #,fill = factor(paste0("Chain_", i))
        geom_line(aes(y = !!sym(paste0("starting_", i))), color = fills_plot[i], linetype = 2, linewidth = 0.7) + #, color = factor(paste0("Chain_", i))
        # Add a line for the 50th percentile (median) for each chain.
        geom_line(aes(y = !!sym(paste0("Chain", i, "_50"))), color = fills_plot[i], linewidth = 1.15) #, color = factor(paste0("Chain_", i))
    }
    # Add additional lines for 'data' and 'starting' values.
    graph_to_plot <- graph_to_plot + 
      geom_line(aes(y = data), color = "green", linetype = 2, linewidth = 1.5) +
      # Set labels and theme options.
      ylab("Beta") + theme(axis.title.x = element_blank(), legend.position = "none") +
      # Manually set color scales for fills and lines.
      coord_cartesian(ylim = c(0, 0.7))
    # print(graph_to_plot)
  } else {
    # The code for non-synthetic data is similar, but with different plot elements.
    graph_to_plot <- ggplot(beta_plot2, aes(x = day))
    for(i in 1:length(valid_chains)){
      graph_to_plot <- graph_to_plot + geom_ribbon(aes(ymin = !!sym(paste0("Chain", i, "_5")),
                                                      ymax = !!sym(paste0("Chain", i, "_95"))),
              fill = fills_plot[i], alpha = 0.5) +
        geom_line(aes(y = !!sym(paste0("starting_", i))), color = fills_plot[i], linetype = 2, linewidth = 0.7) +
        geom_line(aes(y = !!sym(paste0("Chain", i, "_50"))), color = fills_plot[i], linewidth = 1.15)
    }
    
    # Add a line for 'starting' values and vertical lines for specific dates.
    graph_to_plot <- graph_to_plot + 
      geom_vline(xintercept = as.Date(measures$Date), linetype = "dotdash",
                color = measures$Colour, linewidth = 0.5) +
      # Set axis limits, labels, and theme options.
      coord_cartesian(ylim = c(0, 0.6)) +
      ylab("Transmission Rate") + xlab("Date") + theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
                                                      axis.text = element_text(size = 20),
                                                      axis.title = element_text(size = 27), legend.position = "none") +
      ggtitle(paste0("Transmission rate: ", plotTitle)) +
      scale_x_date(breaks = "month", date_labels = "%m/%y") +
      scale_fill_manual(values = fills_plot) + scale_color_manual(values = fills_plot)
  }

  # Check if plotting is enabled ('do_plot' flag).
  if(do_plot){
    # Print the generated plot.
    print(graph_to_plot)
  }
  pdf(file = paste0(file_directory,"plots/", file_name, "_beta.pdf"), paper = "a4r", width = 0, height = 0)
  print(graph_to_plot)
  dev.off()

  t1 <- Sys.time()
  print(paste0("Finished Beta computation: ", t1 - t0))

  remove(beta)
  remove(beta_total)
  remove(beta_plot2)
  remove(provisional)
  remove(graph_to_plot)
  gc()

  print("Starting R0 computation")
  t0 <- Sys.time()
  # Initialize an empty list for R0 values.
  R0 <- list()

  # Conditional calculations based on the type of model.
  if(model_type == "SIR"){
    # Loop over the number of chains.
    for (j in 1:length(valid_chains)) { # 1:chains
      # Calculate R0 for each chain in the SIR model.
      # This involves exponentiating the product of 'time_basis' and a subset of 'accepted_proposals',
      # and then dividing by the first column of 'accepted_proposals'(the gamma parameter).
      R0[[j]] <- exp(time_basis %*% t(accepted_proposals[[j]][, 6:n_var])) / accepted_proposals[[j]][, 1]
    }
  } else {
    # Similar calculations for models other than SIR.
    for (j in 1:length(valid_chains)) { #1:chains
      # Here the calculation does not use double transposition as in the SIR model case.
      R0[[j]] <- exp(time_basis %*% t(accepted_proposals[[j]][, 8:n_var])) / accepted_proposals[[j]][, 2]
    }
  }

  # Initialize a variable for the first chain.
  j <- 1 # valid_chains[1]

  # Calculate quantiles and mean for the first chain and store it in a data frame.
  # This involves applying the quantile function to each row of R0 for the first chain.
  provisional <- cbind(as.data.frame(t(apply(R0[[j]], 1, quantile, probs = c(0.01, 0.5, 0.99), na.rm = TRUE))), 
                      rowMeans(R0[[j]], na.rm = TRUE), exp(time_basis %*% initial.splines[[j]]) / initial.params[[j]]$gamma)
  # Set column names for this data frame.
  colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"), paste0("starting_", j))

  # Initialize 'R0_plot2' (for plotting R0) with the values from 'provisional'.
  R0_plot2 <- provisional

  # Process data for remaining chains.
  for (j in 2:length(valid_chains)){ # valid_chains[-1]
    provisional <- cbind(as.data.frame(t(apply(R0[[j]], 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))), 
                        rowMeans(R0[[j]], na.rm = TRUE), exp(time_basis %*% initial.splines[[j]]) / initial.params[[j]]$gamma)
    colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"), paste0("starting_", j))
    R0_plot2 <- cbind(R0_plot2, provisional)
  }

  # Conditional processing based on whether the data is synthetic or not.
  if(syn_data){
    # If the data is synthetic, further processing is based on the model type.
    if(model_type == "SIR"){
      # For the SIR model, create a data frame with calculated 'data' and 'starting' values,
      # and bind it to the left side of 'R0_plot2'.
      # 'data' (the true generating R0) is calculated by exponentiating a product and then dividing by 0.1,
      # 'starting' (the initial proposal for the MCMC) is similarly calculated and
      # divided by the first element of 'initial.params'.
      R0_plot2 <- cbind(data.frame(
        day = 1:ndata, 
        data = exp(time_basis_0 %*% c(-1.869874, -1.301393, -0.2422027, -1.51097, -3.304498,
                                      -3.09169, -1.568307, -1.570512, -3.447937, -4.521353,
                                      -3.334759, -2.809118)) / 0.1
      ), R0_plot2)
    } else {
      # For other model types, the calculation is similar but uses the second element of 'initial.params'.
      R0_plot2 <- cbind(data.frame(
        day = 1:ndata, 
        data = exp(time_basis_0 %*% c(-1.869874, -1.301393, -0.2422027, -1.51097, -3.304498,
                                      -3.09169, -1.568307, -1.570512, -3.447937, -4.521353,
                                      -3.334759, -2.809118)) / 0.1
      ), R0_plot2)
    }
  } else {
    # If the data is not synthetic, the processing again depends on the model type.
    if(model_type == "SIR"){
      # For the SIR model, bind a data frame with 'day' and 'starting' values to 'R0_plot2'.
      # 'starting' is calculated as before but divided by the first element of 'initial.params'.
      R0_plot2 <- cbind(data.frame(
        day = Incidence_Data$Date 
      ), R0_plot2)
    } else {
      # For other model types, the calculation is similar but uses the second element of 'initial.params'.
      R0_plot2 <- cbind(data.frame(
        day = Incidence_Data$Date 
      ), R0_plot2)
    }
  }

  # Combine all the R0 values from different chains into a single data structure, regardless of data type.
  # This is done by starting with the first chain and then binding the rest of the chains one by one.
  R0_total <- R0[[1]]
  for(iii in 2:valid_chains){ # 2:chains
    R0_total <- cbind(R0_total, R0[[iii]])
  }

  # Conditional plot creation based on whether the data is synthetic.
  if(syn_data){
    # Initialize the ggplot object for synthetic data.
    graph_to_plot_2 <- ggplot(R0_plot2, aes(x = day))
    
    # Loop over the chains to add elements to the plot for each chain.
    for(i in 1:length(valid_chains)){ # valid_chains
      # Add a ribbon (shaded area) for the 5th to 95th percentile range.
      graph_to_plot_2 <- graph_to_plot_2 + geom_ribbon(aes(ymin = !!sym(paste0("Chain", i, "_5")),
                                                      ymax = !!sym(paste0("Chain", i, "_95"))), fill = fills_plot[i], group = paste0("Chain_", i), alpha = 0.5) +
        # Add a line for the median (50th percentile).
        geom_line(aes(y = !!sym(paste0("starting_", i))), color = fills_plot[i], linetype = 2, linewidth = 0.7) +
        geom_line(aes(y = !!sym(paste0("Chain", i, "_50"))), color = fills_plot[i], linewidth = 1.15)
    }
    
    # Add more elements to the plot.
    graph_to_plot_2 <- graph_to_plot_2 + 
      geom_line(aes(y = data), color = "green", linetype = 2, linewidth = 1.5) +
      ylab("R0") + theme(axis.title.x = element_blank(), legend.position = "none") +
      coord_cartesian(ylim = c(0, 7)) #+ 
      #scale_color_manual(values = fills_plot) + scale_fill_manual(values = fills_plot)
  } else {
    # Initialize the ggplot object for non-synthetic data.
    graph_to_plot_2 <- ggplot(R0_plot2, aes(x = day))
    
    # Similar loop for adding elements for each chain.
    for(i in 1:length(valid_chains)){
      graph_to_plot_2 <- graph_to_plot_2 + geom_ribbon(aes(ymin = !!sym(paste0("Chain", i, "_5")),
                                                          ymax = !!sym(paste0("Chain", i, "_95"))), fill  = fills_plot[i], alpha = 0.5) +
        geom_line(aes(y = !!sym(paste0("starting_", i))), color = fills_plot[i], linetype = 2, linewidth = 0.7) +
        geom_line(aes(y = !!sym(paste0("Chain", i, "_50"))), color = fills_plot[i], linewidth = 1.15)
    }
    
    # Add more elements to the plot, including horizontal and vertical lines.
    graph_to_plot_2 <- graph_to_plot_2 + 
      geom_hline(yintercept = 1, linetype = "dotted", color = "black", linewidth = 0.5) +
      geom_vline(xintercept = as.Date(measures$Date), linetype = "dotdash",
                color = measures$Colour, linewidth = 0.5) + 
      coord_cartesian(ylim = c(0, 11)) +
      ylab("Basic Reproduction Number") + xlab("Date") +
      theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 27), legend.position = "none") +
      labs(title = paste0("Basic reproduction number: ", plotTitle)) +
      scale_fill_manual(values = fills_plot) +
      scale_x_date(breaks = "month", date_labels = "%m/%y") +
      scale_color_manual(values = fills_plot)
  }

  # Check if plotting is enabled and print the plot to the R console.
  if(do_plot){
    print(graph_to_plot_2)
  }
  pdf(file = paste0(file_directory,"plots/",file_name, "_R0.pdf"), paper = "a4r", width = 0, height = 0)
  print(graph_to_plot_2)
  dev.off()  

  pdf(file = paste0(file_name, "_R0.pdf"), paper = "a4r", width = 0, height = 0)
  print(graph_to_plot_2)
  dev.off()  

  t1 <- Sys.time()
  print(paste0("Finished R0 computation: ", t1 - t0))

  remove(R0)
  remove(R0_total)
  remove(R0_plot2)
  remove(provisional)
  remove(graph_to_plot_2)
  gc()

  print("Starting Densities computations")

  if(model_type == "SEIR"){
    # Open a PDF device for the first density plot.
    pdf(file = paste0(file_directory, "plots/", file_name, "_density_alpha.pdf"), paper = "a4", width = 0, height = 0)
    
    # Initialize the chain index.
    j <- 1 # valid_chains[1]
    
    # Set parameters for the plot based on whether the data is synthetic.
    if (syn_data){
      # Use predefined mean, standard deviation, and axis limits for synthetic data.
      alpha_mean <- 0.5
      alpha_sd <- 0.05
      x_0 <- 0.35
      x_1 <- 0.65
      y_0 <- 0
      y_1 <- 15
    } else {
      # Use model parameters for real data.
      alpha_mean <- model.parameters$alpha_mean
      alpha_sd <- model.parameters$alpha_sd
      x_0 <- qnorm(0.005, alpha_mean, alpha_sd)
      x_1 <- qnorm(0.995, alpha_mean, alpha_sd)
      y_0 <- 0
      y_1 <- 1.5 * dnorm(qnorm(0.5, alpha_mean, alpha_sd), alpha_mean, alpha_sd)
    }
    
    # Plot the density of the first chain's 'alpha' parameter.
    plot(density(accepted_proposals[[j]][, 1]),  main = substitute(bold(plotTitleVar: alpha), list(plotTitleVar = plotTitle)), xlab = "", ylab = "",
        col = fills_plot[j], xlim = c(x_0, x_1), ylim = c(y_0, y_1),
        cex.axis = 2.15, cex.main = 3)
    
    # Add density lines for other chains.
    for (j in 2:length(valid_chains)){ # valid_chains[-1]
      lines(density(accepted_proposals[[j]][, 1]), col = fills_plot[j])
    }
    
    # Add a theoretical normal distribution line.
    lines(seq(x_0, x_1, length.out = 500), dnorm(seq(x_0, x_1, length.out = 500), alpha_mean, alpha_sd),
          col = "brown4", lwd = 2, lty = 2)
    
    # If synthetic data, add a line at alpha = 0.5.
    if(syn_data){
      lines(c(0.5, 0.5), c(0, 40), col = "green", lwd = 2, lty = 2)
    }
    
    # Add a line for the initial parameter value.
    # Add a line for the initial parameter value.
    for (i in 1:length(valid_chains)){
      lines(c(initial.params[[i]]$alpha, initial.params[[i]]$alpha), c(0, 2 * y_1), col = fills_plot[i], lwd = 0.7, lty = 2)
    }
    
    # Close the first PDF device.
    dev.off()
    
  }

  if (model_type == "SEIR"){
    i_adjustment <- 0
  } else if (model_type == "SIR"){
    i_adjustment <- 1
  }
  # 
  if (gamma.restrictions[1]==0) {
    # Open a PDF device for the first density plot.
    pdf(file = paste0(file_directory, "plots/", file_name, "_density_gamma.pdf"), paper = "a4", width = 0, height = 0)
    
    # Initialize the chain index.
    j <- 1 # valid_chains[1]
    
    # Set parameters for the plot based on whether the data is synthetic.
    if (syn_data){
      # Use predefined mean, standard deviation, and axis limits for synthetic data.
      gamma_mean <- 0.1
      gamma_sd <- 0.01
      x_0 <- 0
      x_1 <- 0.3
      y_0 <- 0
      y_1 <- 45
    } else {
      # Use model parameters for real data.
      alpha_mean <- model.parameters$alpha_mean
      alpha_sd <- model.parameters$alpha_sd
      x_0 <- qnorm(0.005, alpha_mean, alpha_sd)
      x_1 <- qnorm(0.995, alpha_mean, alpha_sd)
      y_0 <- 0
      y_1 <- 1.5 * dnorm(qnorm(0.5, alpha_mean, alpha_sd), alpha_mean, alpha_sd)
    }
    
    # Plot the density of the first chain's 'gamma' parameter.
    plot(density(accepted_proposals[[j]][, 2 - i_adjustment]),  main = substitute(bold(plotTitleVar: gamma), list(plotTitleVar = plotTitle)), xlab = "", ylab = "",
        col = fills_plot[j], xlim = c(x_0, x_1), ylim = c(y_0, y_1),
        cex.axis = 2.15, cex.main = 3)
    
    # Add density lines for other chains.
    for (j in 2:length(valid_chains)){ # valid_chains[-1]
      lines(density(accepted_proposals[[j]][, 2 - i_adjustment]), col = fills_plot[j])
    }
    
    # Add a theoretical normal distribution line.
    lines(seq(x_0, x_1, length.out = 500), dnorm(seq(x_0, x_1, length.out = 500), gamma_mean, gamma_sd),
          col = "brown4", lwd = 2, lty = 2)
    
    # If synthetic data, add a line at alpha = 0.5.
    if(syn_data){
      lines(c(0.1, 0.1), c(0, 40), col = "green", lwd = 2, lty = 2)
    }
    
    # Add a line for the initial parameter value.
    # Add a line for the initial parameter value.
    for (i in 1:length(valid_chains)){
      lines(c(initial.params[[i]]$gamma, initial.params[[i]]$gamma), c(0, 2 * y_1), col = fills_plot[i], lwd = 0.7, lty = 2)
    }
    
    # Close the first PDF device.
    dev.off()
  }

  # Open a PDF device for the third density plot.
  pdf(file = paste0(file_directory, "plots/", file_name, "_density_S0.pdf"), paper = "a4", width = 0, height = 0)

  # Initialize the chain index.
  j <- 1 # valid_chains[1]

  # Set parameters for the plot based on whether the data is synthetic.
  if (syn_data){
    # Use predefined mean, standard deviation, and axis limits for synthetic data.
    # I0_mean <- 10
    # I0_sd <- 1
    # x_0 <- 5
    # x_1 <- 15
    y_0 <- 0
    y_1 <- 0.05
  } else {
    # Use model parameters for real data.
    I0_mean <- model.parameters$I0_mean
    I0_sd <- model.parameters$I0_sd
    x_0 <- qnorm(0.005, I0_mean, I0_sd)
    x_1 <- qnorm(0.995, I0_mean, I0_sd)
    y_0 <- 0
    y_1 <- 1.5 * dnorm(qnorm(0.5, I0_mean, I0_sd), I0_mean, I0_sd)
  }

  # Plot the density of the first chain's 'S0' parameter.
  plot(density(accepted_proposals[[j]][, 3 - i_adjustment]),  main = substitute(bold(plotTitleVar: S[0]), list(plotTitleVar = plotTitle)), xlab = "", ylab = "",
      col = fills_plot[j],  ylim = c(y_0, y_1), # xlim = c(x_0, x_1), 
      cex.axis = 2.15, cex.main = 3)

  # Add density lines for other chains.
  for(j in 2:length(valid_chains)){ # valid_chains[-1]
    lines(density(accepted_proposals[[j]][, 3 - i_adjustment]), col = fills_plot[j])
  }
  # Add a line for the initial parameter value.
  for (i in 1:length(valid_chains)){
    lines(c(initial.params[[i]]$S0, initial.params[[i]]$S0), c(0, 2 * y_1), col = fills_plot[i], lwd = 0.7, lty = 2)
  }

  # # Define the scaled density function
  # scaled_dbeta <- function(y, a, b, N) {
  #   if (y < 0 || y > N) return(0)  # Density is 0 outside [0, N]
  #   N * dbeta(y / N, a, b)  # Scale the density
  # }

  y_values <- npop * rbeta(10000, dirich_alpha[1], sum(dirich_alpha) - dirich_alpha[1])
  lines(density(y_values),
        col = "brown4", lwd = 2, lty = 2)

  # If synthetic data, add a line at 'I0' = 10.
  if(syn_data){
    lines(c(npop - 10, npop - 10), c(0, 40), col = "green", lwd = 2, lty = 2)
  }

  # Close the PDF device.
  dev.off()

  if (model_type == "SEIR"){
    # Open a PDF device for the third density plot.
    pdf(file = paste0(file_directory, "plots/", file_name, "_density_E0.pdf"), paper = "a4", width = 0, height = 0)
    # Set parameters for the plot based on whether the data is synthetic.
    if (syn_data){
      # Use predefined mean, standard deviation, and axis limits for synthetic data.
      # I0_mean <- 10
      # I0_sd <- 1
      # x_0 <- 5
      # x_1 <- 15
      y_0 <- 0
      y_1 <- 0.065
    } else {
      # Use model parameters for real data.
      I0_mean <- model.parameters$I0_mean
      I0_sd <- model.parameters$I0_sd
      x_0 <- qnorm(0.005, I0_mean, I0_sd)
      x_1 <- qnorm(0.995, I0_mean, I0_sd)
      y_0 <- 0
      y_1 <- 1.5 * dnorm(qnorm(0.5, I0_mean, I0_sd), I0_mean, I0_sd)
    }
    j <- 1 # valid_chains[1]
    # Plot the density of the first chain's 'E0' parameter.
    plot(density(accepted_proposals[[j]][, 4 - i_adjustment]), main = substitute(bold(plotTitleVar: E[0]), list(plotTitleVar = plotTitle)), xlab = "", ylab = "",
        col = fills_plot[j], ylim = c(y_0, y_1), # xlim = c(x_0, x_1), 
        cex.axis = 2.15, cex.main = 3)
    
    # Add density lines for other chains.
    for(j in 2:length(valid_chains)){ # valid_chains[-1]
      lines(density(accepted_proposals[[j]][, 4 - i_adjustment]), col = fills_plot[j])
    }
    
    # Add a line for the initial parameter value.
    for (i in 1:length(valid_chains)){ # valid_chains
      lines(c(initial.params[[i]]$E0, initial.params[[i]]$E0), c(0, 2 * y_1), col = fills_plot[i], lwd = 0.7, lty = 2)
    }
    
    y_values <- npop * rbeta(10000, dirich_alpha[2], sum(dirich_alpha) - dirich_alpha[2])
    lines(density(y_values),
          col = "brown4", lwd = 2, lty = 2)
    
    # If synthetic data, add a line at alpha = 0.5.
    if(syn_data){
      lines(c(10, 10), c(0, 40), col = "green", lwd = 2, lty = 2)
    }
    dev.off()
  }
  if (model_type == "SEIR"){
    i_adjustmen <- 0
  } else if (model_type == "SIR"){
    i_adjustment <- 2
  }

  # Open a PDF device for the third density plot.
  pdf(file = paste0(file_directory, "plots/", file_name, "_density_I0.pdf"), paper = "a4", width = 0, height = 0)
  # Set parameters for the plot based on whether the data is synthetic.
  if (syn_data){
    # Use predefined mean, standard deviation, and axis limits for synthetic data.
    # I0_mean <- 10
    # I0_sd <- 1
    # x_0 <- 5
    # x_1 <- 15
    y_0 <- 0
    y_1 <- 0.35
  } else {
    # Use model parameters for real data.
    I0_mean <- model.parameters$I0_mean
    I0_sd <- model.parameters$I0_sd
    x_0 <- qnorm(0.005, I0_mean, I0_sd)
    x_1 <- qnorm(0.995, I0_mean, I0_sd)
    y_0 <- 0
    y_1 <- 1.5 * dnorm(qnorm(0.5, I0_mean, I0_sd), I0_mean, I0_sd)
  }
  j <- 1 # valid_chains[1]
  # Plot the density of the first chain's 'I0' parameter.
  plot(density(accepted_proposals[[j]][, 5 - i_adjustment]), main = substitute(bold(plotTitleVar: I[0]), list(plotTitleVar = plotTitle)), xlab = "", ylab = "",
      col = fills_plot[j], ylim = c(y_0, y_1), # xlim = c(x_0, x_1), 
      cex.axis = 2.15, cex.main = 3)

  # Add density lines for other chains.
  for(j in 2:length(valid_chains)){ # valid_chains[-1]
    lines(density(accepted_proposals[[j]][, 5 - i_adjustment]), col = fills_plot[j])
  }
  # Add a line for the initial parameter value.
  for (i in 1:length(valid_chains)){
    lines(c(initial.params[[i]]$I0, initial.params[[i]]$I0), c(0, 2 * y_1), col = fills_plot[i], lwd = 0.7, lty = 2)
  }
  if (i_adjustment == 0){
    y_values <- npop * rbeta(10000, dirich_alpha[3], sum(dirich_alpha) - dirich_alpha[3])
  } else{
    y_values <- npop * rbeta(10000, dirich_alpha[2], sum(dirich_alpha) - dirich_alpha[2])
  }
  lines(density(y_values),
        col = "brown4", lwd = 2, lty = 2)
  # If synthetic data, add a line at alpha = 0.5.
  if(syn_data){
    lines(c(0, 0), c(0, 40), col = "green", lwd = 2, lty = 2)
  }
  dev.off()

  # Open a PDF device for the third density plot.
  pdf(file = paste0(file_directory, "plots/", file_name, "_density_phi_inv.pdf"), paper = "a4", width = 0, height = 0)

  # Initialize the chain index.
  j <- 1 # valid_chains[1]

  # Set parameters for the plot based on whether the data is synthetic.
  if (syn_data){
    # Use predefined values for synthetic data.
    phi_inv_mean <- 10
    x_0 <- 0
    x_1 <- 0.3
    y_0 <- 0
    y_1 <- 25
  } else {
    # Use model parameters for real data.
    phi_inv_mean <- model.parameters$phi_inv_lambda
    x_0 <- qexp(0.005, phi_inv_mean)
    x_1 <- qexp(0.995, phi_inv_mean)
    y_0 <- 0
    y_1 <- 20 * dexp(1/phi_inv_mean, phi_inv_mean)
  }

  # Plot the density of the first chain's phi_inv parameter.
  plot(density(accepted_proposals[[j]][, 6 - i_adjustment]), main = substitute(bold(plotTitleVar: Phi^-1), list(plotTitleVar = plotTitle)), xlab = "", ylab = "",
      col = fills_plot[j], ylim = c(y_0, y_1), xlim = c(x_0, x_1), # , xlim = c(0.05, 2)
      cex.axis = 2.15, cex.main = 3)

  # Add density lines for other chains.
  for (j in 2:length(valid_chains)) { # valid_chains[-1]
    lines(density(accepted_proposals[[j]][, 6 - i_adjustment]), col = j)
  }

  # Add a theoretical exponential distribution curve.
  curve(dexp(x, phi_inv_mean), xlim = c(x_0, x_1),
        col = "brown4", lwd = 2, lty = 2, add = TRUE)

  # If synthetic data, add a specific line.
  if(syn_data){
    lines(c(0.1, 0.1), c(0, 100), col = "green", lwd = 2, lty = 2)
  }

  # Add a line for the initial parameter value.
  for (i in 1:length(valid_chains)){
    lines(c(initial.params[[i]]$phi_inv, initial.params[[i]]$phi_inv), c(0, 2 * y_1), col = fills_plot[i], lwd = 0.7, lty = 2)
  }
  # Close the PDF device.
  dev.off()
  # 

  # Open a PDF device for the third density plot.
  pdf(file = paste0(file_directory, "plots/", file_name, "_density_tau.pdf"), paper = "a4", width = 0, height = 0)

  # Initialize the chain index.
  j <- 1 # valid_chains[1]

  # Set parameters for the plot based on whether the data is synthetic.
  if (syn_data){
    # Use predefined values for synthetic data.
    phi_inv_mean <- 10
    x_0 <- 0
    x_1 <- 25
    y_0 <- 0
    y_1 <- 0.35
  } else {
    # Use model parameters for real data.
    phi_inv_mean <- model.parameters$phi_inv_lambda
    x_0 <- qexp(0.005, phi_inv_mean)
    x_1 <- qexp(0.995, phi_inv_mean)
    y_0 <- 0
    y_1 <- 20 * dexp(1/phi_inv_mean, phi_inv_mean)
  }

  # Plot the density of the first chain's tau parameter.
  plot(density(accepted_proposals[[j]][, 7 - i_adjustment]), main = substitute(bold(plotTitleVar: tau), list(plotTitleVar = plotTitle)), xlab = "", ylab = "",
      col = fills_plot[j], xlim = c(x_0, x_1), ylim = c(y_0, y_1),
      cex.axis = 2.15, cex.main = 3)

  # Add density lines for other chains.
  for (j in 2:length(valid_chains)) { # valid_chains[-1]
    lines(density(accepted_proposals[[j]][, 7- i_adjustment]), col = fills_plot[j])
  }

  # Add a theoretical inverse gamma distribution curve.
  lines(seq(0, 100, length.out = 500), dinvgamma(seq(0, 100, length.out = 500), shape = 1, rate = 0.005),
        col = "brown4", lwd = 2, lty = 2)

  # Add a line for the initial parameter value.
  for (i in 1:length(valid_chains)){
    lines(c(initial.params[[i]]$tau, initial.params[[i]]$tau), c(0, 2 * y_1), col = fills_plot[i], lwd = 0.7, lty = 2)
  }
  dev.off()
  gc()
} 