# -----------------------------
# Save the initialization results for HAICS
# -----------------------------

saving_initialization_to_haics <- function(inits, N, test_name, input_file_path) {

    for (i in 1:10){
        file_name <- paste0(test_name, "_", i)
        full_dir  <- file.path("/haics/simulation/output/", file_name)
        
        # Create directory if it doesn't exist
        if (!dir.exists(full_dir)) {
            dir.create(full_dir, recursive = TRUE)
        }

        # Write initialpoint
        inits_haics_1 <- c(inits[[i]]$alpha_free, inits[[i]]$gamma_free, (inits[[i]]$theta_free * N_haics)[1:3], 
                            inits[[i]]$phi_inv_free, inits[[i]]$tau)
        cat(paste(inits_haics_1, collapse = " "), " ",
            file = file.path(full_dir, paste0("initialpoint_", file_name, ".txt")), sep = "")

        # Write splinebasis
        inits_haics_2 <- as.numeric(inits[[i]]$beta_free)
        cat(paste(inits_haics_2, collapse = " "), " ",
            file = file.path(full_dir, paste0("splinebasis_", file_name, ".txt")), sep = "")

        # Copy input file and rename it
        new_input_file <- file.path(full_dir, paste0("inputfile_", test_name, "_", i, ".txt"))
        file.copy(from = input_file_path, to = new_input_file, overwrite = TRUE)
    }
}
