AnalysisHaicsResults <- function(file_directory, file_name, chains, burnin, time_basis, n_var,
                                 model_type, neq, npop, ndata, gamma_bounds = c(1/30, 1),
                                 Incidence_data, measures, do_plot = FALSE, syn_data = FALSE,
                                 plotTitle = "", chain_analysis = FALSE){
  # Epidemiological Model Analysis Using MCMC Methods
  # This script is part of HaiCS a project that uses Hamiltonian based Markov Chain Monte Carlo (HMC) methods
  # for statistical and epidemiological analysis,focusing on SIKR/SEMIKR models with a time-dependent
  # transmission rate represented in a B-spline basis.
  
  # ---- Posterior Predictive Checks and Data Processing ----
  # The script performs posterior predictive checks by simulating new data
  # based on the MCMC results and processes these results for further analysis.
  # It calculates quantiles and means for different variables, preparing data for visualization.
  
  # ---- Handling Synthetic and Real Data ----
  # Distinguishes between synthetic and real data in its processing.
  # This distinction is key for model validation (with synthetic data) and real-world data analysis.
  
  # ---- Effective Sample Size (ESS) Calculation ----
  # Calculates the Effective Sample Size (ESS) for each chain.
  # ESS is a critical measure of sample quality in MCMC simulations.
  
  # ---- Diagnostic Plot Creation ----
  # Uses 'ggplot2' to create diagnostic plots like trace plots and Gelman plots.
  # These plots are essential for assessing MCMC chain convergence and behavior.
  
  # ---- Compilation of Key Statistics ----
  # Compiles important statistics from the model, such as (time-dependent) transmission rates and reproduction numbers.
  # This step is crucial for interpreting the model's implications.
  
  # ---- Memory Management ----
  # Includes calls to garbage collection (`gc()`) for efficient memory use.
  # Important in the context of computationally intensive MCMC simulations.
  
  # Note: The script is structured to handle different scenarios and parameters,
  # adapting to either synthetic or real datasets.
  
  # ---- Arguments ----
  # file_directory: A string with the directory where the outputs of HaiCS can be found.
  #                 For a simulation called TEST it will be ".../haics/output/TEST".
  # file_name: A string indicating the name of the simulation ran by HaiCS.
  #            Following the previous example it wold be "TEST".
  # chains: The number of chains that were run with HaiCS.
  # burnin: The amount of burnin that was used in the Hamiltonian-based MCMC sampling.
  # time_basis: a matrix where each column represents a basis function, and each row
  #             corresponds to a time point. The elements of this matrix are the values
  #             of the B-spline basis functions evaluated at these time points.
  # n_var: Number of parameters for which one wants the posterior distribution.
  #        The number of parameters sampled in the Hamiltonian-based MCMC. 
  # neq: Number of equations for the epidemiological compartmental model.
  #      For a SIR model it will be 4 (we add an extra incidence counting equation).
  # npop: The fixed size of the total susceptible population.
  # ndata: The sample size of observed incidence data.
  #       In a daily counting, the number of days for which incidence data is available.
  # gamma_bounds: Upper and lower limits for the gamma (inverse average infectious time) parameter of the
  #               SIKR/SEMIKR model.
  # Incidence_data: A data frame with variables Date (yyyy-mm-dd for the daily incidence data),
  #                Cases (daily new positive cases), Modif_Cases (incidence data corrected for underreporting)
  # measures: A data frame with variables Measure (String indicating the (non-)pharmacological measure taken),
  #          Date (yyyy-mm-dd of the measure), Colour (for plotting pourposes).
  # do_plot: A boolean stating if plots should be shown on screen.
  # syn_data: A boolean stating if synthetic data, for which true parameters are known, have been used.
  #           This triggers plotting the true values for reference.
  # plotTitle: A string meant to indicate the type of compartmental model used.
  # chain_analysis: A boolean indicating if convergence analysis of the chains should be performed.
  #                 If TRUE cooda's gelman and trace plots are performed.
  
  # ---- Value ----
  # A list with the following components:
  # ESS: A matrix with the ESS of each chain (rows) and variable (columns).
  # AR: The acceptance rate for each chain.
  # TransmissionRate: Extracted transmission rate statistics.
  # ReproductionNumber: Extracted reproduction number statistics.
  # PosteriroPredictiveMean: Posterior predictive means.
  # beta: Full beta (time-dependent transmission rate) data used for plotting.
  # R0: Full R0 (time-dependent basic reproduction number) data used for plotting .
  # post_pred: Full posterior predictive data used for plotting.
  # beta_full: Combined beta data across all chains.
  # R0_full: Combined R0 data across all chains.
  # post_pred_full: Combined posterior predictive data across all chains.
  # PDF and SVG files with plots of the transmission rate, basic reproduction number and
  # posterior densities of the epidemiological parameters will be saved in file_directory.
  # Optionally trace and Gelman could be provided.
  
  require(outbreaks)
  require(ggplot2)
  require(splines)
  require(coda)
  require(invgamma)
  require(RColorBrewer)
  require(parallel)
  require(data.table)
  require(MASS)
  
  # Defining the bounds for the gamma parameter 
  gammaL <- gamma_bounds[1]
  gammaU <- gamma_bounds[2]
  
  # Determine the number of cores to use (leaving one free for the system)
  no_cores <- min(max(detectCores() - 1, 1), chains)
  
  # Parallel processing of chains using mclapply
  results <- mclapply(1:chains, function(i) {
    traj_data <- fread(paste0(file_directory,"/",file_name,"_", i,"/trajectories", ".txt"))
    ode_data <- fread(paste0(file_directory,"/",file_name,"_", i,"/ode_sol", ".txt"))
    
    trajectories_original <- if(model_type == "SIR"){
      data.frame(
        Gamma = gammaL + (gammaU - gammaL) / (1 + exp(-traj_data[, 1])),
        I0 = 1 + (npop - 1) / (1 + exp(-traj_data[, 2])),
        Phi_Inv = exp(traj_data[, 3]) + 1e-10,
        Tau = exp(traj_data[, 4]),
        traj_data[, 5:n_var]
      )
    } else {
      data.frame(
        Alpha = exp(traj_data[, 1]) + 1e-10,
        Gamma = gammaL + (gammaU - gammaL) / (1 + exp(-traj_data[, 2])),
        I0 = 1 + (npop - 1) / (1 + exp(-traj_data[, 3])),
        Phi_Inv = exp(traj_data[, 4]) + 1e-10,
        Tau = exp(traj_data[, 5]),
        traj_data[, 6:n_var]
      )
    }
    
    list(trajectories_original = trajectories_original, ode_sol = ode_data)
  }, mc.cores = no_cores)
  
  # Unpack results
  trajectories_original <- lapply(results, `[[`, "trajectories_original")
  ode_sol <- lapply(results, `[[`, "ode_sol")
  
  
  # Create empty lists for accepted proposals and ODE solutions
  accepted_proposals <- list()
  accepted_ode_sol <- list()
  
  # Process each chain for accepted proposals
  for (i in 1:chains) {
    burn_in = burnin
    start_idx = ifelse(burn_in == 0, 1, burn_in + 1)
    
    # Extract the relevant part of trajectories_original after burn-in
    relevant_traj = trajectories_original[[i]][start_idx:nrow(trajectories_original[[i]]), ]
    
    # Get unique row indices from the relevant part of trajectories_original
    unique_indices = which(!duplicated(relevant_traj))
    
    # Adjust indices to align with original ode_sol (since unique_indices are based on post-burn-in subset)
    adjusted_indices = unique_indices + burn_in - 1
    
    # Use the adjusted unique indices to extract rows from trajectories_original and ode_sol
    accepted_proposals[[i]] <- relevant_traj[unique_indices, ]
    accepted_ode_sol[[i]] <- ode_sol[[i]][adjusted_indices, ]
  }

  # Create an array to store the dimensions of accepted proposals for each chain
  dim_a_prop <- array(0, dim = chains)
  
  # Create an array to store the acceptance rate for each chain
  acceptance_rate <- array(0, dim = chains)
  
  # Loop over each chain
  for (i in 1:chains) {
    # Calculate the number of rows (dimensions) of accepted proposals for the current chain
    dim_a_prop[i] <- dim(accepted_proposals[[i]])[1]
    
    # Calculate the acceptance rate for the current chain
    acceptance_rate[i] <- dim_a_prop[i] / (dim(trajectories_original[[i]])[1] - burnin)
  }
  
  # We assume that initial splines and parameters are the same for each chain
  # for clarity in the plots.
  
  # Read and store the data from a file containing initial splines
  initial.splines <- t(read.table(paste0(file_directory, "/input/splinebasis_", file_name, "_1.txt"), quote="\"", comment.char=""))
  
  # Read and store the data from a file containing initial parameters
  initial.params <- t(read.table(paste0(file_directory, "/input/initialpoint_", file_name, "_1.txt"), quote="\"", comment.char=""))
  
  # Check if the number of chains is 10 or less.
  if (chains <= 10){
    # If there are 10 or fewer chains, set a color palette using the 'brewer.pal' function.
    # The palette is set to have 10 colors from the "Paired" palette.
    palette(brewer.pal(n = 10, name = "Paired"))
  } else {
    # If there are more than 10 chains, set a color palette with a number of colors
    # equal to the number of chains. This ensures each chain gets a unique color.
    palette(brewer.pal(n = chains, name = "Paired"))
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
  
  # Initialize an empty list called 'beta' (the time-dependent transmission rate).
  beta <- list()
  
  # Check the type of model being used.
  if(model_type == "SIR"){
    # If the model type is "SIR", perform a calculation for each chain.
    for (j in 1:chains) {
      # Calculate 'beta' for each chain using a matrix product ('%*%') between 'time_basis'
      # and a subset of 'accepted_proposals', then exponentiate the result.
      beta[[j]] <- exp(time_basis %*% t(accepted_proposals[[j]][, 5:n_var]))
    }
  } else {
    # If the model type is not "SIR", use a similar calculation but start from the 6th column.
    for (j in 1:chains) {
      beta[[j]] <- exp(time_basis %*% t(accepted_proposals[[j]][, 6:n_var]))
    }
  }
  
  # Initialize a variable for the first chain.
  j <- 1
  
  # Calculate quantiles and mean for the first chain and store it in a data frame.
  # The 'apply' function is used to calculate quantiles and mean for each row in 'beta[[j]]'.
  provisional <- cbind(as.data.frame(t(apply(beta[[j]], 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))), 
                       rowMeans(beta[[j]], na.rm = TRUE))
  # Set column names for the 'provisional' dataframe.
  colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"))
  
  # Initialize 'beta_plot2' (which will be used for plotting the transmissionrate)
  # with the values from 'provisional'.
  beta_plot2 <- provisional
  
  # Process data for remaining chains.
  for (j in 2:chains){
    provisional <- cbind(as.data.frame(t(apply(beta[[j]], 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))), 
                         rowMeans(beta[[j]], na.rm = TRUE))
    colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"))
    beta_plot2 <- cbind(beta_plot2, provisional)
  }
  
  # Conditional plotting based on whether the data is synthetic ('syn_data') or not.
  if (syn_data){
    # If the data is synthetic, bind additional columns (the true transmission rate) to 'beta_plot2'.
    beta_plot2 <- cbind(data.frame(
      day = 1:ndata, 
      data = exp(time_basis %*% c(-1.869874, -1.301393, -0.2422027, -1.51097, -3.304498, -3.09169, -1.568307, -1.570512, -3.447937, -4.521353, -3.334759, -2.809118)),
      starting = exp(time_basis %*% initial.splines)
    ), beta_plot2)
  } else {
    # If the data is not synthetic, bind different columns (the initial proposal for the MCMC) to 'beta_plot2'.
    beta_plot2 <- cbind(data.frame(
      day = Incidence_data$Date, 
      starting = exp(time_basis %*% initial.splines)
    ), beta_plot2)
  }
  
  # Combine all chains into one data structure.
  beta_total <- beta[[1]]
  for(iii in 2:chains){ 
    beta_total <- cbind(beta_total, beta[[iii]])
  }

  # Conditional plot generation based on whether the data is synthetic or not.
  if(syn_data){
    # Start building a plot with 'day' on the x-axis.
    graph_to_plot <- ggplot(beta_plot2, aes(x = day))
    
    # Loop over the number of chains to add graphical elements for each chain.
    for(i in 1:chains){
      # Add a ribbon (shaded area) for the 5th and 95th percentiles.
      graph_to_plot <- graph_to_plot + geom_ribbon(aes(ymin = get(paste0("Chain", i, "_5")),
                                                       ymax = get(paste0("Chain", i, "_95")), 
                                                       fill = paste0("Chain_", i)), alpha = 0.5) +
        # Add a line for the 50th percentile (median) for each chain.
        geom_line(aes(y = get(paste0("Chain", i, "_50")), color = paste0("Chain_", i)), linewidth = 1.15)
    }
    
    # Add additional lines for 'data' and 'starting' values.
    graph_to_plot <- graph_to_plot + 
      geom_line(aes(y = data), color = "green", linetype = 2, linewidth = 1.5) +
      geom_line(aes(y = starting), color = "gray30", linetype = 2, linewidth = 1.5) +
      # Set labels and theme options.
      ylab("Beta") + theme(axis.title.x = element_blank(), legend.position = "none") +
      # Manually set color scales for fills and lines.
      scale_fill_manual(values = fills_plot) + scale_color_manual(values = fills_plot)
  } else {
    # The code for non-synthetic data is similar, but with different plot elements.
    graph_to_plot <- ggplot(beta_plot2, aes(x = day))
    for(i in 1:chains){
      graph_to_plot <- graph_to_plot + geom_ribbon(aes(ymin = get(paste0("Chain", i, "_5")),
                                                       ymax = get(paste0("Chain", i, "_95")), 
                                                       fill = paste0("Chain_", i)), alpha = 0.5) +
        geom_line(aes(y = get(paste0("Chain", i, "_50")), color = paste0("Chain_", i)), linewidth = 1.15)
    }
    
    # Add a line for 'starting' values and vertical lines for specific dates.
    graph_to_plot <- graph_to_plot + 
      geom_line(aes(y = starting), color = "gray30", linetype = 2, linewidth = 0.7) +
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
  
  # Save the plot as a PDF file.
  # 'file_directory' and 'file_name' are variables that store the path to the directory and the base name for the file.
  # 'paste0' is used to concatenate these strings along with "_beta.pdf" to create the full file name.
  # 'paper = "a4r"' specifies the paper size as A4 in landscape orientation (‘a4r’ for A4 rotated).
  # 'width = 0' and 'height = 0' imply that the default size of the PDF device will be used.
  pdf(file = paste0(file_directory, "/", file_name, "_beta.pdf"), paper = "a4r", width = 0, height = 0)
  
  # Print the plot to the PDF device. This is the step where the plot gets drawn into the file.
  print(graph_to_plot)
  
  # Close the PDF graphics device. This step is necessary to complete the writing process to the file.
  dev.off()
  
  # Save the plot as an SVG file.
  # Similar to the PDF, the path and file name are constructed using 'paste0'.
  # The dimensions of the SVG are explicitly set here: 7.8 inches in height and 11.2 inches in width.
  svg(file = paste0(file_directory, "/", file_name, "_beta.svg"), height = 7.8, width = 11.2)
  
  # Print the plot to the SVG device.
  print(graph_to_plot)
  
  # Close the SVG graphics device.
  dev.off()
  
  # Initialize an empty list for R0 values.
  R0 <- list()
  
  # Conditional calculations based on the type of model.
  if(model_type == "SIR"){
    # Loop over the number of chains.
    for (j in 1:chains) {
      # Calculate R0 for each chain in the SIR model.
      # This involves exponentiating the product of 'time_basis' and a subset of 'accepted_proposals',
      # and then dividing by the first column of 'accepted_proposals'(the gamma parameter).
      R0[[j]] <- exp(time_basis %*% t(accepted_proposals[[j]][, 5:n_var])) / accepted_proposals[[j]][, 1]
    }
  } else {
    # Similar calculations for models other than SIR.
    for (j in 1:chains) {
      # Here the calculation does not use double transposition as in the SIR model case.
      R0[[j]] <- exp(time_basis %*% t(accepted_proposals[[j]][, 6:n_var])) / accepted_proposals[[j]][, 2]
    }
  }
  
  # Initialize a variable for the first chain.
  j <- 1
  
  # Calculate quantiles and mean for the first chain and store it in a data frame.
  # This involves applying the quantile function to each row of R0 for the first chain.
  provisional <- cbind(as.data.frame(t(apply(R0[[j]], 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))), 
                       rowMeans(R0[[j]], na.rm = TRUE))
  # Set column names for this data frame.
  colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"))
  
  # Initialize 'R0_plot2' (for plotting R0) with the values from 'provisional'.
  R0_plot2 <- provisional
  
  # Process data for remaining chains.
  for (j in 2:chains){
    provisional <- cbind(as.data.frame(t(apply(R0[[j]], 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))), 
                         rowMeans(R0[[j]], na.rm = TRUE))
    colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"))
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
        data = exp(time_basis %*% c(-1.869874, -1.301393, -0.2422027, -1.51097, -3.304498,
                                    -3.09169, -1.568307, -1.570512, -3.447937, -4.521353,
                                    -3.334759, -2.809118)) / 0.1, 
        starting = exp(time_basis %*% initial.splines) / initial.params[1]
      ), R0_plot2)
    } else {
      # For other model types, the calculation is similar but uses the second element of 'initial.params'.
      R0_plot2 <- cbind(data.frame(
        day = 1:ndata, 
        data = exp(time_basis %*% c(-1.869874, -1.301393, -0.2422027, -1.51097, -3.304498,
                                    -3.09169, -1.568307, -1.570512, -3.447937, -4.521353,
                                    -3.334759, -2.809118)) / 0.1, 
        starting = exp(time_basis %*% initial.splines) / initial.params[2]
      ), R0_plot2)
    }
  } else {
    # If the data is not synthetic, the processing again depends on the model type.
    if(model_type == "SIR"){
      # For the SIR model, bind a data frame with 'day' and 'starting' values to 'R0_plot2'.
      # 'starting' is calculated as before but divided by the first element of 'initial.params'.
      R0_plot2 <- cbind(data.frame(
        day = Incidence_data$Date, 
        starting = exp(time_basis %*% initial.splines) / initial.params[1]
      ), R0_plot2)
    } else {
      # For other model types, the calculation is similar but uses the second element of 'initial.params'.
      R0_plot2 <- cbind(data.frame(
        day = Incidence_data$Date, 
        starting = exp(time_basis %*% initial.splines) / initial.params[2]
      ), R0_plot2)
    }
  }
  
  # Combine all the R0 values from different chains into a single data structure, regardless of data type.
  # This is done by starting with the first chain and then binding the rest of the chains one by one.
  R0_total <- R0[[1]]
  for(iii in 2:chains){
    R0_total <- cbind(R0_total, R0[[iii]])
  }
  
  # Conditional plot creation based on whether the data is synthetic.
  if(syn_data){
    # Initialize the ggplot object for synthetic data.
    graph_to_plot <- ggplot(R0_plot2, aes(x = day))
    
    # Loop over the chains to add elements to the plot for each chain.
    for(i in 1:chains){
      # Add a ribbon (shaded area) for the 5th to 95th percentile range.
      graph_to_plot <- graph_to_plot + geom_ribbon(aes(ymin = get(paste0("Chain", i, "_5")),
                                                       ymax = get(paste0("Chain", i, "_95")), fill = paste0("Chain_", i)), alpha = 0.5) +
        # Add a line for the median (50th percentile).
        geom_line(aes(y = get(paste0("Chain", i, "_50")), color = paste0("Chain_", i)), linewidth = 1.15)
    }
    
    # Add more elements to the plot.
    graph_to_plot <- graph_to_plot + 
      geom_line(aes(y = data), color = "green", linetype = 2, linewidth = 1.5) +
      geom_line(aes(y = starting), color = "gray30", linetype = 2, linewidth = 0.7) +
      ylab("R0") + theme(axis.title.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = fills_plot) + coord_cartesian(ylim = c(0, 45)) + 
      scale_color_manual(values = fills_plot)
  } else {
    # Initialize the ggplot object for non-synthetic data.
    graph_to_plot_2 <- ggplot(R0_plot2, aes(x = day))
    
    # Similar loop for adding elements for each chain.
    for(i in 1:chains){
      graph_to_plot_2 <- graph_to_plot_2 + geom_ribbon(aes(ymin = get(paste0("Chain", i, "_5")),
                                                           ymax = get(paste0("Chain", i, "_95")), fill = paste0("Chain_", i)), alpha = 0.5) +
        geom_line(aes(y = get(paste0("Chain", i, "_50")), color = paste0("Chain_", i)), linewidth = 1.15)
    }
    
    # Add more elements to the plot, including horizontal and vertical lines.
    graph_to_plot_2 <- graph_to_plot_2 + 
      geom_line(aes(y = starting), color = "gray30", linetype = 2, linewidth = 0.7) +
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
  
  # Save the plot to a PDF file.
  pdf(file = paste0(file_directory, "/", file_name, "_R0.pdf"), paper = "a4r", width = 0, height = 0)
  # Print the plot to the PDF device.
  print(graph_to_plot_2)
  # Close the PDF graphics device.
  dev.off()
  
  # Check if plotting is enabled and print the plot to the R console.
  if(do_plot){
    print(graph_to_plot_2)
  }
  
  # Save the plot to a SVG file.
  svg(file = paste0(file_directory,"/",file_name, "_R0.svg"), height = 7.8, width = 11.2)
  # Print the plot to the SVG device.
  print(graph_to_plot_2)
  # Close the PDF graphics device.
  dev.off()
  
  # Plotting posterior parameter densities conditionally on the model
  if(model_type == "SEIR"){
    # Read model parameters from a text file. We are assuming all chains have the same input_file.
    model.parameters.0 <- read.table(paste0(file_directory, "/input/inputfile_", file_name, "_1.txt"),
                                     quote="\"", comment.char="")
    
    # Process the model parameters into a data frame.
    suppressWarnings(model.parameters <- data.frame(t(as.numeric(model.parameters.0[, 2]))))
    names(model.parameters) <- model.parameters.0[, 1]
    
    # Open a PDF device for the first density plot.
    pdf(file = paste0(file_directory, "/", file_name, "_density_1.pdf"), paper = "a4", width = 0, height = 0)
    
    # Initialize the chain index.
    j = 1
    
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
    plot(density(accepted_proposals[[j]][, 1]), main = expression(bold(paste(plotTitle, ": ", alpha))), xlab = "", ylab = "",
         col = j, xlim = c(x_0, x_1), ylim = c(y_0, y_1),
         cex.axis = 2.15, cex.main = 3)
    
    # Add density lines for other chains.
    for (j in 2:chains){
      lines(density(accepted_proposals[[j]][, 1]), col = j)
    }
    
    # Add a theoretical normal distribution line.
    lines(seq(x_0, x_1, length.out = 500), dnorm(seq(x_0, x_1, length.out = 500), alpha_mean, alpha_sd),
          col = "brown4", lwd = 2, lty = 2)
    
    # If synthetic data, add a line at alpha = 0.5.
    if(syn_data){
      lines(c(0.5, 0.5), c(0, 40), col = "green", lwd = 2, lty = 2)
    }
    
    # Add a line for the initial parameter value.
    lines(c(initial.params[1], initial.params[1]), c(0, 2 * y_1), col = "gray30", lwd = 0.7, lty = 2)
    
    # Close the first PDF device.
    dev.off()
    # Open a PDF device for the second density plot.
    pdf(file = paste0(file_directory, "/", file_name, "_density_2.pdf"), paper = "a4", width = 0, height = 0)
    
    # Initialize the chain index.
    j = 1
    
    # Set parameters for the plot based on whether the data is synthetic.
    if (syn_data){
      # Use predefined mean, standard deviation, and axis limits for synthetic data.
      gamma_mean <- 0.1
      gamma_sd <- 0.05
      x_0 <- 0
      x_1 <- 0.20
      y_0 <- 0
      y_1 <- 20
    } else {
      # Use model parameters for real data.
      gamma_mean <- model.parameters$gamma_mean
      gamma_sd <- model.parameters$gamma_sd
      x_0 <- 1/30  # A predefined lower limit (the reciprocal of an assumed maximum recovery time).
      x_1 <- qnorm(0.995, gamma_mean, gamma_sd)  # Upper limit based on the normal distribution.
      y_0 <- 0
      y_1 <- 1.5 * dnorm(qnorm(0.5, gamma_mean, gamma_sd), gamma_mean, gamma_sd)
    }
    
    # Plot the density of the first chain's 'gamma' parameter.
    plot(density(accepted_proposals[[j]][, 2]), main = expression(bold(paste(plotTitle, ": ", gamma))),
         xlab = "", ylab = "", col = j, xlim = c(x_0, x_1), ylim = c(y_0, y_1),
         cex.axis = 2.15, cex.main = 3)
    
    # Add density lines for other chains.
    for (j in 2:chains){
      lines(density(accepted_proposals[[j]][, 2]), col = j)
    }
    
    # Add a theoretical normal distribution curve.
    curve(dnorm(x, gamma_mean, gamma_sd), xlim = c(x_0, x_1),
          col = "brown4", lwd = 2, lty = 2, add = TRUE)
    
    # Uncomment the following line to add a uniform distribution curve.
    # curve(dunif(x, x_0, x_1),
    #       col = "brown4", lwd = 2, lty = 2, add = TRUE)
    
    # If synthetic data, add a line at gamma = 0.1.
    if(syn_data){
      lines(c(0.1, 0.1), c(0, 40), col = "green", lwd = 2, lty = 2)
    }
    
    # Add a line for the initial parameter value.
    lines(c(initial.params[2], initial.params[2]), c(0, 2 * y_1), col = "gray30", lwd = 0.7, lty = 2)
    
    # Close the PDF device.
    dev.off()
    
    # Open a PDF device for the third density plot.
    pdf(file = paste0(file_directory, "/", file_name, "_density_3.pdf"), paper = "a4", width = 0, height = 0)
    
    # Initialize the chain index.
    j = 1
    
    # Set parameters for the plot based on whether the data is synthetic.
    if (syn_data){
      # Use predefined mean, standard deviation, and axis limits for synthetic data.
      I0_mean <- 10
      I0_sd <- 1
      x_0 <- 5
      x_1 <- 15
      y_0 <- 0
      y_1 <- 0.6
    } else {
      # Use model parameters for real data.
      I0_mean <- model.parameters$I0_mean
      I0_sd <- model.parameters$I0_sd
      x_0 <- qnorm(0.005, I0_mean, I0_sd)
      x_1 <- qnorm(0.995, I0_mean, I0_sd)
      y_0 <- 0
      y_1 <- 1.5 * dnorm(qnorm(0.5, I0_mean, I0_sd), I0_mean, I0_sd)
    }
    
    # Plot the density of the first chain's 'I0' parameter.
    plot(density(accepted_proposals[[j]][, 3]), main = expression(bold(paste(plotTitle, ": ", I[0]))), xlab = "", ylab = "",
         col = j, xlim = c(x_0, x_1), ylim = c(y_0, y_1),
         cex.axis = 2.15, cex.main = 3)
    
    # Add density lines for other chains.
    for(j in 2:chains){
      lines(density(accepted_proposals[[j]][, 3]), col = j)
    }
    
    # Add a theoretical normal distribution curve.
    curve(dnorm(x, I0_mean, I0_sd), xlim = c(x_0, x_1),
          col = "brown4", lwd = 2, lty = 2, add = TRUE)
    
    # If synthetic data, add a line at 'I0' = 10.
    if(syn_data){
      lines(c(10, 10), c(0, 40), col = "green", lwd = 2, lty = 2)
    }
    
    # Add a line for the initial parameter value.
    lines(c(initial.params[3], initial.params[3]), c(0, 2 * y_1), col = "gray30", lwd = 0.7, lty = 2)
    
    # Close the PDF device.
    dev.off()
    
    # Open a PDF device for the fourth density plot (phi_inv parameter).
    pdf(file = paste0(file_directory, "/", file_name, "_density_4.pdf"), paper = "a4", width = 0, height = 0)
    
    # Initialize the chain index.
    j = 1
    
    # Set parameters for the plot based on whether the data is synthetic.
    if (syn_data){
      # Use predefined values for synthetic data.
      phi_inv_mean <- 10
      x_0 <- 0.05
      x_1 <- 0.2
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
    plot(density(accepted_proposals[[j]][, 4]), main = expression(bold(paste(plotTitle, ": ", Phi^-1))), xlab = "", ylab = "",
         col = j, xlim = c(x_0, x_1), ylim = c(y_0, y_1),
         cex.axis = 2.15, cex.main = 3)
    
    # Add density lines for other chains.
    for (j in 2:chains) {
      lines(density(accepted_proposals[[j]][, 4]), col = j)
    }
    
    # Add a theoretical exponential distribution curve.
    curve(dexp(x, phi_inv_mean), xlim = c(x_0, x_1),
          col = "brown4", lwd = 2, lty = 2, add = TRUE)
    
    # If synthetic data, add a specific line.
    if(syn_data){
      lines(c(0.1, 0.1), c(0, 100), col = "green", lwd = 2, lty = 2)
    }
    
    # Add a line for the initial parameter value.
    lines(c(initial.params[4], initial.params[4]), c(0, 2 * y_1), col = "gray30", lwd = 0.7, lty = 2)
    
    # Close the PDF device.
    dev.off()
    
    # Open a PDF device for the fifth density plot (tau parameter).
    pdf(file = paste0(file_directory, "/", file_name, "_density_5.pdf"), paper = "a4", width = 0, height = 0)
    
    # Initialize the chain index.
    j = 1
    
    # Plot the density of the first chain's tau parameter.
    plot(density(accepted_proposals[[j]][, 5]), main = expression(bold(paste(plotTitle, ": ", tau))), xlab = "", ylab = "",
         col = j, xlim = c(0, 100), ylim = c(0, 0.08),
         cex.axis = 2.15, cex.main = 3)
    
    # Add density lines for other chains.
    for (j in 2:chains) {
      lines(density(accepted_proposals[[j]][, 5]), col = j)
    }
    
    # Add a theoretical inverse gamma distribution curve.
    lines(seq(0, 100, length.out = 500), dinvgamma(seq(0, 100, length.out = 500), shape = 1, rate = 0.005),
          col = "brown4", lwd = 2, lty = 2)
    
    # Add a line for the initial parameter value.
    lines(c(initial.params[4], initial.params[4]), c(0, 100), col = "gray30", lwd = 2, lty = 2)
    
    # Close the PDF device.
    dev.off()
  }
  else{ # The same as for the SEIR model but for a SIR model
    suppressWarnings(model.parameters <- data.frame(t(as.numeric(model.parameters.0[, 2]))))
    names(model.parameters) <- model.parameters.0[, 1]
    
    pdf(file = paste0(file_directory,"/",file_name, "_density_1.pdf"), paper = "a4", width = 0, height = 0)
    j = 1
    if (syn_data){
      gamma_mean <- 0.1
      gamma_sd <- 0.05
      x_0 <- 0
      x_1 <- 0.20
      y_0 <- 0
      y_1 <- 20
    } else{
      gamma_mean <- model.parameters$gamma_mean
      gamma_sd <- model.parameters$gamma_sd
      x_0 <- 1/30
      x_1 <- qnorm(0.995, model.parameters$gamma_mean, model.parameters$gamma_sd)
      y_0 <- 0
      y_1 <- 1.5 * dnorm(qnorm(0.5, model.parameters$gamma_mean, model.parameters$gamma_sd),
                         model.parameters$gamma_mean, model.parameters$gamma_sd)
    }
    plot(density(accepted_proposals[[j]][, 1]), main = expression(bold(paste(plotTitle,": ", gamma))),
         xlab = "", ylab = "", col = j, xlim = c(x_0, x_1), ylim = c(y_0, y_1),
         cex.axis = 2.15, cex.main = 3)
    for (j in 2:chains){
      lines(density(accepted_proposals[[j]][, 1]), col = j)
    }
    curve(dnorm(x, gamma_mean, gamma_sd), xlim = c(x_0, x_1),
          col = "brown4", lwd = 2, lty = 2, add = TRUE)
    # curve(dunif(x, x_0, x_1),
    #       col = "brown4", lwd = 2, lty = 2, add = TRUE)
    if(syn_data){
      lines(c(0.1, 0.1), c(0, 40),
            col = "green", lwd = 2, lty = 2)
    }
    lines(c(initial.params[1], initial.params[1]),  c(0, 2 * y_1),
          col = "gray30", lwd = 0.7, lty = 2)
    dev.off()
    pdf(file = paste0(file_directory,"/",file_name, "_density_2.pdf"), paper = "a4", width = 0, height = 0)
    j = 1
    if (syn_data){
      I0_mean <- 10
      I0_sd <- 1
      x_0 <- 5
      x_1 <- 15
      y_0 <- 0
      y_1 <- 0.4
    } else{
      I0_mean <- model.parameters$I0_mean
      I0_sd <- model.parameters$I0_sd
      x_0 <- qnorm(0.005, model.parameters$I0_mean, model.parameters$I0_sd)
      x_1 <- qnorm(0.995, model.parameters$I0_mean, model.parameters$I0_sd)
      y_0 <- 0
      y_1 <- 1.5 * dnorm(qnorm(0.5, model.parameters$I0_mean, model.parameters$I0_sd),
                         model.parameters$I0_mean, model.parameters$I0_sd)
    }
    plot(density(accepted_proposals[[j]][, 2]), main = expression(bold(paste(plotTitle,": ", I[0]))),
         xlab = "", ylab = "",
         col = j, xlim = c(x_0, x_1), ylim = c(y_0, y_1),
         cex.axis = 2.15, cex.main = 3)
    for(j in 2:chains){
      lines(density(accepted_proposals[[j]][, 2]), col = j)
    }
    curve(dnorm(x, I0_mean, I0_sd), xlim = c(x_0, x_1),
          col = "brown4", lwd = 2, lty = 2, add = TRUE)
    if(syn_data){
      lines(c(10, 10), c(0, 40),
            col = "green", lwd = 2, lty = 2)
    }
    lines(c(initial.params[2], initial.params[2]),  c(0, 2 * y_1),
          col = "gray30", lwd = 0.7, lty = 2)
    dev.off()
    pdf(file = paste0(file_directory,"/",file_name, "_density_3.pdf"), paper = "a4", width = 0, height = 0)
    j = 1
    if (syn_data){
      phi_inv_mean <- 10
      x_0 <- 0.05
      x_1 <- 0.2
      y_0 <- 0
      y_1 <- 25
    } else{
      phi_inv_mean <- model.parameters$phi_inv_lambda
      x_0 <- qexp(0.005, model.parameters$phi_inv_lambda)
      x_1 <- qexp(0.995, model.parameters$phi_inv_lambda)
      y_0 <- 0
      y_1 <- 20 * dexp(1/model.parameters$phi_inv_lambda, model.parameters$phi_inv_lambda)
    }
    plot(density(accepted_proposals[[j]][, 3]), main = expression(bold(paste(plotTitle,": ", Phi^-1))),
         xlab = "", ylab = "", col = j, xlim = c(x_0, x_1), ylim = c(y_0, y_1),
         cex.axis = 2.15, cex.main = 3)
    for (j in 2:chains) {
      lines(density(accepted_proposals[[j]][, 3]), col = j)
    }
    curve(dexp(x, phi_inv_mean), xlim = c(x_0, x_1),
          col = "brown4", lwd = 2, lty = 2, add = TRUE)
    if(syn_data){
      lines(c(0.1, 0.1), c(0, 100),
            col = "green", lwd = 2, lty = 2)
    }
    lines( c(initial.params[3], initial.params[3]),  c(0, 2 * y_1),
           col = "gray30", lwd = 0.7, lty = 2)
    dev.off()
    
    pdf(file = paste0(file_directory,"/",file_name, "_density_4.pdf"), paper = "a4", width = 0, height = 0)
    j = 1
    plot(density(accepted_proposals[[j]][, 4]), main = expression(bold(paste(plotTitle,": ", tau))),
         xlab = "", ylab = "", col = j, xlim = c(0, 100), ylim = c(0, 0.15),
         cex.axis = 2.15, cex.main = 3)
    for (j in 2:chains) {
      lines(density(accepted_proposals[[j]][, 4]), col = j)
    }
    lines(seq(0, 100, length.out = 500), dinvgamma(seq(0, 100, length.out = 500), shape = 1, rate = 0.005),
          col = "brown4", lwd = 2, lty = 2)
    lines(c(initial.params[4], initial.params[4]),  c(0, 100),
          col = "gray30", lwd = 2, lty = 2)
    dev.off()
  }
  
  # Check if chain analysis is to be performed.
  if(chain_analysis){
    # Check if plotting is enabled.
    if (do_plot){
      # Set graphical parameters based on whether the data is synthetic.
      if (syn_data){
        par(mfrow = c(6, 3), mai = c(0.26, 0.26, 0.2, 0.2))
      } else{
        par(mfrow = c(7, 3), mai = c(0.26, 0.26, 0.2, 0.2))
      }
      
      # Create an empty list to store the MCMC objects.
      mcmc_objects <- list()
      
      # Convert trajectories data into MCMC objects for each chain.
      for (i in 1:chains) {
        mcmc_objects[[i]] <- as.mcmc(trajectories_original[[i]][-(1:burn_in),])
      }
      
      # Plot Gelman diagnostics.
      gelman.plot(mcmc.list(mcmc_objects), auto.layout = FALSE)
      
      # Reset graphical parameters for trace plots.
      if (syn_data){
        par(mfrow = c(6, 3), mai = c(0.26, 0.26, 0.2, 0.2))
      } else{
        par(mfrow = c(7, 3), mai = c(0.26, 0.26, 0.2, 0.2))
      }
      
      # Generate trace plots for the MCMC objects.
      coda::traceplot(mcmc_objects, col = palette())
    }
    
    # Save Gelman diagnostics plot to a PDF file.
    pdf(file = paste0(file_directory, "/", file_name, "_gelman.pdf"), paper = "a4", width = 0, height = 0)
    # Set graphical parameters.
    if (syn_data){
      par(mfrow = c(6, 3), mai = c(0.26, 0.26, 0.2, 0.2))
    } else{
      par(mfrow = c(7, 3), mai = c(0.26, 0.26, 0.2, 0.2))
    }
    # Generate and save the Gelman plot.
    gelman.plot(mcmc.list(mcmc_objects), auto.layout = FALSE)
    # Close the PDF device for the Gelman plot.
    dev.off()
    
    # Save trace plots to a PDF file.
    pdf(file = paste0(file_directory, "/", file_name, "_trace.pdf"), paper = "a4", width = 0, height = 0)
    # Set graphical parameters.
    if (syn_data){
      par(mfrow = c(6, 3), mai = c(0.26, 0.26, 0.2, 0.2))
    } else{
      par(mfrow = c(7, 3), mai = c(0.26, 0.26, 0.2, 0.2))
    }
    # Generate and save the trace plot.
    coda::traceplot(mcmc_objects, col = palette())
    # Close the PDF device for the trace plot.
    dev.off()
  }

  # Initialize an array to store the Effective Sample Size (ESS) for each chain and each variable.
  # 'chains' is the number of MCMC chains, and 'n_var' is the number of variables.
  ESS_chains <- array(dim = c(chains, n_var))
  
  # Loop over each chain.
  for (i in 1:chains){
    # Convert the trajectory of the i-th chain to an MCMC object, excluding the burn-in period.
    # 'burn_in' is the number of initial steps of the chain to be discarded.
    # 'trajectories_original[[i]]' contains the trajectory of the i-th chain.
    mcmc_obj <- as.mcmc(trajectories_original[[i]][(burn_in+1):dim(trajectories_original[[i]])[1], ])
    
    # Calculate the Effective Sample Size for each variable in the i-th chain.
    # The 'effectiveSize' function is part of the 'coda' package and calculates ESS.
    ESS_chains[i, ] <- effectiveSize(mcmc_obj)
  }
  
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
  posterior_predictive_check <- mclapply(1:chains, function(j) {
    posterior <- array(dim = c(dim_a_prop[[j]], ndata))
    
    for (i in 1:dim_a_prop[[j]]) {
      # Determine the 'size' parameter based on the model type
      size_param = if(model_type == "SIR") trajectories_original[[j]][i, 3] else trajectories_original[[j]][i, 4]
      
      # Vectorized generation of random numbers using rnegbin
      posterior[i, ] <- rnegbin(as.numeric(model_incidence[[j]][i, ]), theta = 1 / size_param)
    }
    return(posterior)
  }, mc.cores = no_cores)

  # Initialize a data frame for the first chain's posterior predictive check results.
  j <- 1
  provisional <- cbind(
    # Calculate the 5th, 50th (median), and 95th percentiles and the mean for each time point.
    as.data.frame(t(apply(posterior_predictive_check[[j]], 2, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE))),
    colMeans(posterior_predictive_check[[j]], na.rm = TRUE)
  )
  colnames(provisional) <- c(paste0("Chain", j, "_5"), paste0("Chain", j, "_50"), paste0("Chain", j, "_95"), paste0("Chain", j, "_Mean"))
  posterior_predictive_check_plot2 <- provisional
  
  # Repeat the process for the remaining chains.
  for (j in 2:chains){
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
        day = Incidence_data$Date, 
        data_corr = Incidence_data$Modif_Cases[1:ndata], 
        data = Incidence_data$Cases[1:ndata]
      ), 
      posterior_predictive_check_plot2
    )
  }
  
  # Initialize the total posterior predictive check data with the first chain's data.
  posterior_predictive_check_total <- posterior_predictive_check[[1]]
  
  # Loop over the remaining chains to append their posterior predictive check data to the total data.
  for(iii in 2:chains){
    posterior_predictive_check_total <- rbind(posterior_predictive_check_total, posterior_predictive_check[[iii]])
  }
  
  # Conditional plotting based on whether the data is synthetic.
  if(syn_data){
    # Initialize the ggplot object for synthetic data.
    graph_to_plot_3 <- ggplot(posterior_predictive_check_plot2, aes(x = day))
    
    # Loop over the chains to add graphical elements for each chain.
    for(i in 1:chains){
      # Add a ribbon for the 5th to 95th percentile range and a line for the median.
      graph_to_plot_3 <- graph_to_plot_3 + 
        geom_ribbon(aes(ymin = get(paste0("Chain", i, "_5")),
                        ymax = get(paste0("Chain", i, "_95")), fill = paste0("Chain_", i)), alpha = 0.5) +
        geom_line(aes(y = get(paste0("Chain", i, "_50")), color = paste0("Chain_", i)), linewidth = 1.15)
    }
    
    # Add more elements to the plot.
    graph_to_plot_3 <- graph_to_plot_3 + 
      geom_point(aes(y = data)) +
      geom_line(aes(y = Mean), color = "green", linetype = 2, linewidth = 1.5) +
      ylab("Daily Incidence") + theme(axis.title.x = element_blank(), legend.position = "none") +
      scale_fill_manual(values = fills_plot) + scale_color_manual(values = fills_plot)
  } else {
    # Initialize the ggplot object for non-synthetic data.
    graph_to_plot_3 <- ggplot(posterior_predictive_check_plot2, aes(x = day))
    
    # Similar loop for adding elements for each chain.
    for(i in 1:chains){
      graph_to_plot_3 <- graph_to_plot_3 + 
        geom_ribbon(aes(ymin = get(paste0("Chain", i, "_5")),
                        ymax = get(paste0("Chain", i, "_95")), fill = paste0("Chain_", i)), alpha = 0.5) +
        geom_line(aes(y = get(paste0("Chain", i, "_50")), color = paste0("Chain_", i)), linewidth = 1.15)
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
  
  # Save the plot as a PDF file.
  pdf(file = paste0(file_directory, "/", file_name, "_posteriorpredictive.pdf"), paper = "a4r", width = 0, height = 0)
  print(graph_to_plot_3)
  dev.off()
  
  # Save the plot as an SVG file.
  svg(file = paste0(file_directory, "/", file_name, "_posteriorpredictive.svg"), height = 7.8, width = 11.2)
  print(graph_to_plot_3)
  dev.off()
  
  # Return a list of various measures and results for real data.
  return(list(
    ESS = ESS_chains,  # Effective Sample Size for each chain.
    AR = acceptance_rate,  # Acceptance rate of the MCMC algorithm.
    TransmissionRate = beta_plot2[, seq(6, (4 * chains + 2), 4)],  # Extracted transmission rate statistics.
    ReproductionNumber = R0_plot2[, seq(6, (4 * chains + 2), 4)],  # Extracted reproduction number statistics.
    PosteriroPredictiveMean = posterior_predictive_check_plot2[, seq(6, (4 * chains + 2), 4)],  # Posterior predictive means.
    beta = beta_plot2,  # Full beta data used for plotting.
    R0 = R0_plot2,  # Full R0 data used for plotting .
    post_pred = posterior_predictive_check_plot2,  # Full posterior predictive data used for plotting.
    beta_full = beta_total,  # Combined beta data across all chains.
    R0_full = R0_total,  # Combined R0 data across all chains.
    post_pred_full = posterior_predictive_check_total  # Combined posterior predictive data across all chains.
  ))
  
  # Call the garbage collector to clean up memory.
  gc()
}
