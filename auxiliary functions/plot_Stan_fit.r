
plot_Stan_fit <- function(fit_file, data) {
  source("auxiliary function/spline_basis_functions.R") # load the bbase() function for splines
  source("auxiliary function/creating_time_basis.R") # load the time basis creation code
  y <- data

  # ──────────────────────────────────────────────────────────────────────────────
  # 1) Define a color palette for later plots (Okabe–Ito is colorblind‐friendly)
  # ──────────────────────────────────────────────────────────────────────────────

  okabe_ito_10 <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#882255", "#44AA99"
  )

  # ──────────────────────────────────────────────────────────────────────────────
  # 2) Pick which Stan fit file to load
  # ──────────────────────────────────────────────────────────────────────────────

  # Read the fitted Stan model object from disk:
  fit <- readRDS(fit_file)

  gamma_fixed <- fit@inits[[1]]$gamma_free  # fixed gamma value for scaling
  T_total <- length(fit@inits[[1]]$y_rep) # total number of days in the model
  M_splines <- length(fit@inits[[1]]$beta_eff) # number of spline basis functions

  # Summarize the posterior draws (means, sds, Rhat, n_eff, quantiles…)
  fit_summary <- rstan::summary(fit)

  # ──────────────────────────────────────────────────────────────────────────────
  # 3) Examine sampler diagnostics: leapfrog steps per iteration
  #    (HMC/NUTS makes multiple “leapfrog” steps in each iteration)
  # ──────────────────────────────────────────────────────────────────────────────

  # Extract the sampler parameters (one matrix per chain).  inc_warmup=FALSE
  # means we drop warmup/adaptation iterations.
  sampler_params_list <- rstan::get_sampler_params(fit, inc_warmup = FALSE)

  # Pull out just the n_leapfrog__ column
  # (number of integrator steps per iteration)
  n_leapfrog_per_chain <- lapply(sampler_params_list, function(chain_mat) {
    chain_mat[, "n_leapfrog__"]
  })

  # Combine all chains into one vector and show a quick summary
  all_steps <- unlist(n_leapfrog_per_chain)
  summary(all_steps)
  # If these vary wildly or hit max_treedepth often,
  # you may need to adjust control.

  # ──────────────────────────────────────────────────────────────────────────────
  # 3.5) Plot the posterior predictive distribution of new daily cases y_rep[t]
  #    using the merged chains summary (median + 95% CI)
  # ──────────────────────────────────────────────────────────────────────────────
  # derive the directory name from the fit_file basename
  out_dir <- substr(fit_file, 1, nchar(fit_file) - 4)

  # if it doesn’t exist, make it (including any parent paths)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # now open the PDF in that folder
  out_file <- file.path(out_dir,
                        "posterior_predictive_daily_incidence_merged_chains")

  summary_y_pred <- fit_summary$summary[paste0("y_rep[", 1:T_total, "]"), ]
  pdf(paste0(out_file, ".pdf"), width = 8, height = 6)
  plot(1:T_total,
      summary_y_pred[, "50%"],  # column named "50%" is the median
      type = "l",
      ylim = c(0, 450),
      xlab = "Day",
      ylab = "Incidence",
      main = paste0("Posterior predictive daily incidence\n(merged chains, gamma = ", gamma_fixed, ")"))

  # Add dashed red lines for the 2.5% and 97.5% quantiles
  # (columns "2.5%" and "97.5%")
  lines(1:T_total, summary_y_pred[, "2.5%"],  col = "red", lty = 2)
  lines(1:T_total, summary_y_pred[, "97.5%"], col = "red", lty = 2)
  # Overlay the actual observed counts as blue dots
  points(1:T_total, y, pch = 16, col = "blue")

  # Legend to explain colors/linetypes
  legend("topright",
        legend = c("Observed", "Posterior median", "95% CI"),
        col    = c("blue",     "black",           "red"),
        pch    = c(16,         NA,                NA),
        lty    = c(NA,         1,                 2),
        bty    = "n")
  dev.off()

  png(paste0(out_file, ".png"),
    width    = 8,       # in inches
    height   = 6,       # in inches
    units    = "in",    # so width/height are inches, not pixels
    res      = 300      # 300 dpi for a crisp figure)
  )
  plot(1:T_total,
      summary_y_pred[, "50%"],  # column named "50%" is the median
      type = "l",
      ylim = c(0, 450),
      xlab = "Day",
      ylab = "Incidence",
      main = paste0("Posterior predictive daily incidence\n(merged chains, gamma = ", gamma_fixed, ")"))

  # Add dashed red lines for the 2.5% and 97.5% quantiles
  # (columns "2.5%" and "97.5%")
  lines(1:T_total, summary_y_pred[, "2.5%"],  col = "red", lty = 2)
  lines(1:T_total, summary_y_pred[, "97.5%"], col = "red", lty = 2)
  # Overlay the actual observed counts as blue dots
  points(1:T_total, y, pch = 16, col = "blue")

  # Legend to explain colors/linetypes
  legend("topright",
        legend = c("Observed", "Posterior median", "95% CI"),
        col    = c("blue",     "black",           "red"),
        pch    = c(16,         NA,                NA),
        lty    = c(NA,         1,                 2),
        bty    = "n")
  dev.off()

  # ──────────────────────────────────────────────────────────────────────────────
  # 4) Plot the posterior predictive distribution of new daily cases C_pred[t]
  #    using the merged chains summary (median + 95% CI)
  # ──────────────────────────────────────────────────────────────────────────────
  # Extract rows "C_pred[1]" through "C_pred[T_total]" from the fit_summary table:
  summary_C_pred <- fit_summary$summary[paste0("C_pred[", 1:T_total, "]"), ] # nolint

  # now open the PDF in that folder
  out_file <- file.path(out_dir,
                        "posterior_average_daily_incidence_merged_chains")

  # Base plot: x=day, y=posterior median (column 6 of summary),
  # with y‐limits chosen to cover your data range.
  pdf(paste0(out_file, ".pdf"), width = 8, height = 6)
  plot(1:T_total,
      summary_C_pred[, "50%"],  # column named "50%" is the median
      type = "l",
      ylim = c(0, 450),
      xlab = "Day",
      ylab = "Incidence",
      main = paste0("Posterior average daily incidence\n(merged chains, gamma = ", gamma_fixed, " )")
  )

  # Add dashed red lines for the 2.5% and 97.5% quantiles
  # (columns "2.5%" and "97.5%")
  lines(1:T_total, summary_C_pred[, "2.5%"],  col = "red", lty = 2)
  lines(1:T_total, summary_C_pred[, "97.5%"], col = "red", lty = 2)
  # Overlay the actual observed counts as blue dots
  points(1:T_total, y, pch = 16, col = "blue")

  # Legend to explain colors/linetypes
  legend("topright",
        legend = c("Observed", "Posterior median", "95% CI"),
        col    = c("blue",     "black",           "red"),
        pch    = c(16,         NA,                NA),
        lty    = c(NA,         1,                 2),
        bty    = "n")
  dev.off()

  png(paste0(out_file, ".png"),
    width    = 8,       # in inches
    height   = 6,       # in inches
    units    = "in",    # so width/height are inches, not pixels
    res      = 300      # 300 dpi for a crisp figure)
  )
  plot(1:T_total,
      summary_C_pred[, "50%"],  # column named "50%" is the median
      type = "l",
      ylim = c(0, 450),
      xlab = "Day",
      ylab = "Incidence",
      main = paste0("Posterior average daily incidence\n(merged chains, gamma = ", gamma_fixed, " )"))

  # Add dashed red lines for the 2.5% and 97.5% quantiles
  # (columns "2.5%" and "97.5%")
  lines(1:T_total, summary_C_pred[, "2.5%"],  col = "red", lty = 2)
  lines(1:T_total, summary_C_pred[, "97.5%"], col = "red", lty = 2)
  # Overlay the actual observed counts as blue dots
  points(1:T_total, y, pch = 16, col = "blue")

  # Legend to explain colors/linetypes
  legend("topright",
        legend = c("Observed", "Posterior median", "95% CI"),
        col    = c("blue",     "black",           "red"),
        pch    = c(16,         NA,                NA),
        lty    = c(NA,         1,                 2),
        bty    = "n")
  dev.off()

  # ──────────────────────────────────────────────────────────────────────────────
  # 5) Plot chain‐wise trajectories of the same C_pred[t]
  #    to check for chain agreement or mixing issues.
  # ──────────────────────────────────────────────────────────────────────────────

  # Re‐extract with chain dimension: this gives an array
  summary_C_pred_chainwise <- rstan::summary( # nolint
    fit,
    pars = paste0("C_pred[", 1:T_total, "]"),
    use_cache = FALSE
  )$c_summary  # this is an array [T_total, quantile, chain]

  # now open the PDF in that folder
  out_file <- file.path(out_dir, "posterior_average_daily_incidence_by_chain")

  pdf(paste0(out_file, ".pdf"), width = 8, height = 6)
  # First, plot the median of chain 1:
  plot(1:T_total,
      summary_C_pred_chainwise[, "50%", 1],
      type = "l",
      ylim = c(0, 450),
      xlab = "Day",
      ylab = "Incidence",
      main = paste0("Posterior pred. by chain (10 chains, gamma = ", gamma_fixed, ")"))

  # Overlay data
  points(1:T_total, y, pch = 16, col = "blue")

  # For each chain, plot median and 95% CI in a different color
  for (chain in 1:10) {
    col_i <- okabe_ito_10[chain]
    med  <- summary_C_pred_chainwise[, "50%", chain]
    lo   <- summary_C_pred_chainwise[, "2.5%", chain]
    hi   <- summary_C_pred_chainwise[, "97.5%", chain]
    lines(1:T_total, med, col = col_i, lty = 1)
    lines(1:T_total, lo,  col = col_i, lty = 2)
    lines(1:T_total, hi,  col = col_i, lty = 2)
  }

  legend("topright",
    legend = c(paste0("Chain ", 1:10),
              "Median",
              "95% CI"),
    col    = c(okabe_ito_10, "black", "black"),
    lty    = c(rep(1, 10),      1,  2),
    lwd    = c(rep(1, 10),      2,  1),
    ncol   = 2,
    bty    = "n",
    inset  = 0.02,
    # ---- new spacing tweaks ----
    y.intersp  = 1.3,                # increase vertical space
    x.intersp  = 0.7,                # give more room around the text
    text.width = max(strwidth(paste0("Chain ", 1:10))),
    # uniform column width
    cex        = 0.85                # slightly smaller text overall
  )
  dev.off()

  png(paste0(out_file, ".png"),
    width    = 8,       # in inches
    height   = 6,       # in inches
    units    = "in",    # so width/height are inches, not pixels
    res      = 300      # 300 dpi for a crisp figure))
  )# First, plot the median of chain 1:
  plot(1:T_total,
      summary_C_pred_chainwise[, "50%", 1],
      type = "l",
      ylim = c(0, 450),
      xlab = "Day",
      ylab = "Incidence",
      main = paste0("Posterior pred. by chain (10 chains, gamma = ", gamma_fixed, ")"))

  # Overlay data
  points(1:T_total, y, pch = 16, col = "blue")

  # For each chain, plot median and 95% CI in a different color
  for (chain in 1:10) {
    col_i <- okabe_ito_10[chain]
    med  <- summary_C_pred_chainwise[, "50%", chain]
    lo   <- summary_C_pred_chainwise[, "2.5%", chain]
    hi   <- summary_C_pred_chainwise[, "97.5%", chain]
    lines(1:T_total, med, col = col_i, lty = 1)
    lines(1:T_total, lo,  col = col_i, lty = 2)
    lines(1:T_total, hi,  col = col_i, lty = 2)
  }

  legend("topright",
    legend = c(paste0("Chain ", 1:10),
              "Median",
              "95% CI"),
    col    = c(okabe_ito_10, "black", "black"),
    lty    = c(rep(1, 10),      1,  2),
    lwd    = c(rep(1, 10),      2,  1),
    ncol   = 2,
    bty    = "n",
    inset  = 0.02,
    # ---- new spacing tweaks ----
    y.intersp  = 1.3,                # increase vertical space
    x.intersp  = 0.7,                # give more room around the text
    text.width = max(strwidth(paste0("Chain ", 1:10))),  # uniform column width
    cex        = 0.85                # slightly smaller text overall
  )
  dev.off()
  # ──────────────────────────────────────────────────────────────────────────────
  # 6) Check convergence diagnostics for key parameters
  #    (Rhat near 1 and n_eff adequate)
  # ──────────────────────────────────────────────────────────────────────────────

  # Print Rhat and effective sample size for S0,E0,I0 and 3 transmission params
  fit_summary$summary[c("S0", "E0", "I0", "R0"),
                      c("Rhat", "n_eff")]

  fit_summary$summary[c("alpha_eff",
                        "gamma_eff", "phi_inv_eff"),
                      c("Rhat", "n_eff")]

  # And for all M_splines spline coefficients beta_eff[1..M_splines]
  fit_summary$summary[paste0("beta_eff[", 1:M_splines, "]"),
                      c("Rhat", "n_eff")]

  write.csv(fit_summary$summary[c("S0", "E0", "I0", "R0",
                                  "alpha_eff", "gamma_eff", "phi_inv_eff",
                                  paste0("beta_eff[", 1:M_splines, "]")),
                      c("Rhat", "n_eff")],
    file = file.path(out_dir, "stan_fit_summary.csv"),
    row.names = TRUE
  )

  # ──────────────────────────────────────────────────────────────────────────────
  # 7) Plot marginal densities and traceplots for parameters of interest
  # ──────────────────────────────────────────────────────────────────────────────
  library(ggplot2)
  # Density plots
  rstan::stan_dens(fit, pars = c("S0", "E0", "I0"))
  rstan::stan_dens(fit, pars = c("alpha_eff", "gamma_eff", "phi_inv_eff"))

  beta_generator <- c(-1.8699, -1.3014,  -0.2422, -1.5110,
                      -3.3045, -3.0917,  -1.5683, -1.5705,
                      -3.4479, -4.5214,  -3.334,  -2.8091)
  # 1. build the base stan_dens plot
  p <- rstan::stan_dens(fit, pars = paste0("beta_eff[", 1:M_splines, "]"))

  # 2. figure out the *name* of the facet‐variable
  facet_var <- names(p$facet$params$facets)
  #> e.g. facet_var == "variable"

  # 3. build your line‐data, naming that column exactly the same
  line_df <- data.frame(x_line = beta_generator)
  line_df[[facet_var]] <- paste0("beta_eff[", 1:M_splines, "]")

  # 4. add a geom_vline layer; because line_df *does* have the facet‐var
  #    the lines will only show up in the matching panels
  # now open the PDF in that folder
  out_file <- file.path(out_dir, "posterior_beta_merged_chains")

  pdf(paste0(out_file, ".pdf"), width = 14, height = 10)
  print(p +
    geom_vline(
      data        = line_df,
      mapping     = aes(xintercept = x_line),
      inherit.aes = FALSE,
      linetype    = "dashed",
      colour      = "blue",
      size        = 1.5               # thicker line
    ) +
    labs(
      x     = expression(beta[eff]),
      title = paste0("Posterior densities with true generator values (gamma = ", gamma_fixed, " )")
    ) +
    theme_minimal(base_size = 16) +   # increase base text size for axes, titles, etc.
    theme(
      plot.title   = element_text(face = "bold", size = 22),
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 18),
      strip.text   = element_text(size = 18),  # facet labels
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 18)
    )
  )
  dev.off()
  png(paste0(out_file, ".png"),
    width    = 14,       # in inches
    height   = 10,       # in inches
    units    = "in",    # so width/height are inches, not pixels
    res      = 300      # 300 dpi for a crisp figure))
  )
  print(p +
    geom_vline(
      data        = line_df,
      mapping     = aes(xintercept = x_line),
      inherit.aes = FALSE,
      linetype    = "dashed",
      colour      = "blue",
      size        = 1.5               # thicker line
    ) +
    labs(
      x     = expression(beta[eff]),
      title = paste0("Posterior densities with true generator values (gamma = ", gamma_fixed, " )")
    ) +
    theme_minimal(base_size = 16) +   # increase base text size for axes, titles, etc.
    theme(
      plot.title   = element_text(face = "bold", size = 22),
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 18),
      strip.text   = element_text(size = 18),  # facet labels
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 18)
    )
  )
  dev.off()

  out_file <- file.path(out_dir, "trceplot_init_comp_by_chain")
  pdf(paste0(out_file, ".pdf"), width = 18, height = 12)
  # Traceplots (overplot all chains, color by chain)
  print(rstan::traceplot(fit, pars = c("S0", "E0", "I0", "R0")) +
    ggplot2::scale_color_manual(values = okabe_ito_10) +
    theme(
      plot.title   = element_text(face = "bold", size = 22),
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 18),
      strip.text   = element_text(size = 18),  # facet labels
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 18)
    ))
  dev.off()
  png(paste0(out_file, ".png"),
    width    = 14,       # in inches
    height   = 10,       # in inches
    units    = "in",    # so width/height are inches, not pixels
    res      = 300      # 300 dpi for a crisp figure))
  )
  # Traceplots (overplot all chains, color by chain)
  print(rstan::traceplot(fit, pars = c("S0", "E0", "I0", "R0")) +
    ggplot2::scale_color_manual(values = okabe_ito_10) +
    theme(
      plot.title   = element_text(face = "bold", size = 22),
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 18),
      strip.text   = element_text(size = 18),  # facet labels
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 18)
    ))
  dev.off()

  rstan::traceplot(fit, pars = c("alpha_eff", "gamma_eff", "phi_inv_eff")) +
    ggplot2::scale_color_manual(values = okabe_ito_10)

  out_file <- file.path(out_dir, "trceplot_beta_by_chain")

  pdf(paste0(out_file, ".pdf"), width = 18, height = 12)
  print(rstan::traceplot(fit,
                  pars = paste0("beta_eff[", 1:M_splines, "]")) +
    ggplot2::scale_color_manual(values = okabe_ito_10) +
    theme(
      plot.title   = element_text(face = "bold", size = 22),
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 18),
      strip.text   = element_text(size = 18),  # facet labels
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 18)
    )
  )
  dev.off()
  png(paste0(out_file, ".png"),
    width    = 14,       # in inches
    height   = 10,       # in inches
    units    = "in",    # so width/height are inches, not pixels
    res      = 300      # 300 dpi for a crisp figure))
  )
  print(rstan::traceplot(fit,
                  pars = paste0("beta_eff[", 1:M_splines, "]")) +
    ggplot2::scale_color_manual(values = okabe_ito_10) +
    theme(
      plot.title   = element_text(face = "bold", size = 22),
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 18),
      strip.text   = element_text(size = 18),  # facet labels
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 18)
    )
  )
  dev.off()

  ; # Autocorrelation plots
  ; rstan::stan_ac(fit, pars = c("S0", "E0", "I0", "R0"))
  ; rstan::stan_ac(fit, pars = c("alpha_eff", "gamma_eff", "phi_inv_eff"))


  # 1. Pull out the lp__ draws: a matrix [iterations × chains]
  lp_array <- as.array(fit)[, , "lp__"]

  out_file <- file.path(out_dir, "log_posterior_traceplot")

  # 2. Traceplot of lp__ for each chain
  pdf(paste0(out_file, ".pdf"), width = 8, height = 6)
  matplot(lp_array, type = "l", lty = 1, col = okabe_ito_10,
          xlab = "Iteration", ylab = "log-posterior",
          main = paste0("Trace of log-posterior by chain (gamma = ", gamma_fixed, ")"))
  legend("topright", legend = paste("chain", seq_len(ncol(lp_array))),
        col = okabe_ito_10, lty = 1, cex = 0.8)
  dev.off()

  png(paste0(out_file, ".png"),
    width    = 8,       # in inches
    height   = 6,       # in inches
    units    = "in",    # so width/height are inches, not pixels
    res      = 300      # 300 dpi for a crisp figure))
  )
  matplot(lp_array, type = "l", lty = 1, col = okabe_ito_10,
          xlab = "Iteration", ylab = "log-posterior",
          main = paste0("Trace of log-posterior by chain (gamma = ", gamma_fixed, ")"))
  legend("topright", legend = paste("chain", seq_len(ncol(lp_array))),
        col = okabe_ito_10, lty = 1, cex = 0.8)
  dev.off()

  # 1. Get the sampler params list (one matrix per chain)
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = TRUE)

  # 2. Pull out the energy__ column for each chain
  energy_list <- lapply(sampler_params, function(ch) ch[, "energy__"])

  # 3. Combine into a matrix [iter × chain]
  energy_mat <- do.call(cbind, energy_list)
  out_file <- file.path(out_dir, "energy_traceplot")
  # 4. Traceplot
  pdf(paste0(out_file, ".pdf"), width = 8, height = 6)
  matplot(energy_mat, type = "l", lty = 1, col = okabe_ito_10,
          xlab = "Iteration", ylab = "HMC energy",
          main = paste0("Trace of HMC energy by chain (gamma = ", gamma_fixed, ")"), ylim = c(-300, -200))
  legend("topright", legend = paste("chain", seq_len(ncol(energy_mat))),
        col = okabe_ito_10, lty = 1, cex = 0.8)
  dev.off()

  png(paste0(out_file, ".png"),
    width    = 8,       # in inches
    height   = 6,       # in inches
    units    = "in",    # so width/height are inches, not pixels
    res      = 300      # 300 dpi for a crisp figure))
  )
  matplot(energy_mat, type = "l", lty = 1, col = okabe_ito_10,
          xlab = "Iteration", ylab = "HMC energy",
          main = paste0("Trace of HMC energy by chain (gamma = ", gamma_fixed, ")"), ylim = c(-300, -200))
  legend("topright", legend = paste("chain", seq_len(ncol(energy_mat))),
        col = okabe_ito_10, lty = 1, cex = 0.8)
  dev.off()


  ##############################
  # Warmup diagnostics
  ##################################

  # how many warmup iterations per chain?
  warmup_iters <- fit@sim$warmup

  # pull back all sampler params (including warmup)
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = TRUE)

  # set up a 2×(nChains/2) plot grid
  n_chains <- length(sampler_params)

  pdf(file = paste0(out_dir, "/sampler_epsilon_traceplots.pdf"),
      width = 8, height = 6)
  par(mfrow = c(2, n_chains / 2), mar = c(3, 3, 2, 1))

  for (chain in 1:10) {
    sp <- sampler_params[[chain]]
    # 1. Step‐size trace
    plot(sp[, "stepsize__"], type = "l",
        main = paste0("Chain ", chain, ": ε"),
        xlab = "Iteration", ylab = "stepsize",
        ylim = c(0, 0.5))
    abline(v = warmup_iters, col = "red", lty = 2)
  }
  dev.off()
  # compare final step‐sizes across chains
  final_eps <- sapply(sampler_params, function(sp) tail(sp[, "stepsize__"], 1))
  print(final_eps)

  times <- rstan::get_elapsed_time(fit) / 3600 # convert seconds to hours
  total_times <- rowSums(times)
  grand_total <- max(total_times)
  print(times)
  print(total_times) # convert seconds to hours
  print(paste0("Total run time: ", round(grand_total, 2), " hours"))


  # install.packages(c("posterior", "dplyr"))
  library(posterior)
  library(dplyr)

  # 1. Get the draws as a 3-D array [iteration, chain, parameter]
  pars       <- paste0("beta_eff[", 1:M_splines, "]")
  draws_arr  <- as.array(fit, pars = pars)
  # dims(draws_arr) = c(iterations, chains, M_splines)

  n_iter     <- dim(draws_arr)[1]
  n_chain    <- dim(draws_arr)[2]
  n_summary  <- T_total   # your function returns T_total values

  # 2. Preallocate a 3-D array for your summaries:
  #    [iteration, chain, summary_index]
  summary_arr <- array(NA_real_,
                      dim = c(n_iter, n_chain, n_summary),
                      dimnames = list(
                        iteration = NULL,
                        chain     = paste0("chain", seq_len(n_chain)),
                        summary_i = paste0("stat", seq_len(n_summary))
                      ))

  # 3. Your custom function, for illustration:
  summary_beta <- function(betas) {
    # betas is length-M_splines numeric; returns T_total-length numeric
    exp(time_basis_0 %*% betas)
  }

  # 4. Fill in summary_arr
  for (ch in seq_len(n_chain)) {
    # apply over iterations: for chain ch, draws_arr[,ch,] is [iter × M_splines]
    # apply row-wise, getting a matrix [T_total × iter]
    tmp <- apply(draws_arr[ , ch, ], 1, summary_beta)
    # transpose to [iter × T_total] and store
    summary_arr[ , ch, ] <- t(tmp)
  }

  # Now df_summary has columns:
  #   iteration, chain, stat_index (1–T_total), stat_value
  # which you can group/summarize however you like.
  pooled_summary <- matrix(
    summary_arr,
    nrow = dim(summary_arr)[1] * dim(summary_arr)[2],
    ncol = dim(summary_arr)[3]
  )

  out_file <- file.path(out_dir, "posterior_beta_merged_chains")

  pdf(paste0(out_file, ".pdf"), width = 12, height = 10)
  plot(1:T_total, apply(pooled_summary, 2, median), typ = "l",
      ylim = c(0, 1),
      xlab = "Days",
      ylab = "Transmission rate",
      main = paste0("Pooled transmission rate across chains (gamma = ", gamma_fixed, " )"))
  lines(1:T_total, apply(pooled_summary, 2, quantile, 0.025), col = "red", lty = 2)
  lines(1:T_total, apply(pooled_summary, 2, quantile, 0.975), col = "red", lty = 2)
  lines(1:T_total, exp(time_basis_0 %*% beta_generator), col = "blue", lwd = 2)
  # Legend to explain colors/linetypes
  legend("topright",
        legend = c("Posterior median", "95% CI", "True generator"),
        col    = c("black", "red", "blue"),
        lty    = c(1, 2, 1),
        bty    = "n")
  dev.off()

  # 1. Get the draws as a 3-D array [iteration, chain, parameter]
  pars       <- c(paste0("beta_eff[", 1:M_splines, "]"), "gamma_eff")
  draws_arr  <- as.array(fit, pars = pars)
  # dims(draws_arr) = c(iterations, chains, M_splines)

  n_iter     <- dim(draws_arr)[1]
  n_chain    <- dim(draws_arr)[2]
  n_summary  <- T_total   # your function returns 100 values

  # 2. Preallocate a 3-D array for your summaries:
  #    [iteration, chain, summary_index]
  summary_arr <- array(NA_real_,
                      dim = c(n_iter, n_chain, n_summary),
                      dimnames = list(
                        iteration = NULL,
                        chain     = paste0("chain", seq_len(n_chain)),
                        summary_i = paste0("stat", seq_len(n_summary))
                      ))

  # 3. Your custom function, for illustration:
  summary_R0 <- function(params) {
    # betas is length-M_splines numeric; returns 100-length numeric
    exp(time_basis_0 %*% params[1:M_splines]) / params[13]
  }

  # 4. Fill in summary_arr
  for (ch in seq_len(n_chain)) {
    # apply over iterations: for chain ch, draws_arr[,ch,] is [iter × M_splines]
    # apply row-wise, getting a matrix [100 × iter]
    tmp <- apply(draws_arr[ , ch, ], 1, summary_R0)
    # transpose to [iter × 100] and store
    summary_arr[ , ch, ] <- t(tmp)
  }

  # now open the PDF in that folder
  out_file <- file.path(out_dir, "posterior_R0_by_chain")

  pdf(paste0(out_file, ".pdf"), width = 8, height = 6)
  # First, plot the median of chain 1:
  col_i <- okabe_ito_10[1]
  med  <- apply(summary_arr[, 1, ], 2, median)
  lo   <- apply(summary_arr[, 1, ], 2, quantile, 0.025)
  hi   <- apply(summary_arr[, 1, ], 2, quantile, 0.975)
  plot(1:T_total, med, col = col_i, type = "l",
      ylim = c(0, 10),
      xlab = "Day",
      ylab = "Basic reproduction number (R0)",
      main = paste0("Posterior R0 by chain (gamma = ", gamma_fixed, ")"), cex = 2)
  lines(1:T_total, lo,  col = col_i, lty = 2)
  lines(1:T_total, hi,  col = col_i, lty = 2)


  # For each chain, plot median and 95% CI in a different color
  for (chain in 2:10) {
    col_i <- okabe_ito_10[chain]
    med  <- apply(summary_arr[, chain, ], 2, median)
    lo   <- apply(summary_arr[, chain, ], 2, quantile, 0.025)
    hi   <- apply(summary_arr[, chain, ], 2, quantile, 0.975)
    lines(1:T_total, med, col = col_i, lty = 1)
    lines(1:T_total, lo,  col = col_i, lty = 2)
    lines(1:T_total, hi,  col = col_i, lty = 2)
  }

  legend("topright",
    legend = c(paste0("Chain ", 1:10),
              "Median",
              "95% CI"),
    col    = c(okabe_ito_10, "black", "black"),
    lty    = c(rep(1, 10),      1,  2),
    lwd    = c(rep(1, 10),      2,  1),
    ncol   = 2,
    bty    = "n",
    inset  = 0.02,
    # ---- new spacing tweaks ----
    y.intersp  = 1.3,                # increase vertical space
    x.intersp  = 0.7,                # give more room around the text
    text.width = max(strwidth(paste0("Chain ", 1:10))),
    # uniform column width
    cex        = 0.85                # slightly smaller text overall
  )
  dev.off()

  # Now df_summary has columns:
  #   iteration, chain, stat_index (1–T_total), stat_value
  # which you can group/summarize however you like.
  pooled_summary <- matrix(
    summary_arr,
    nrow = dim(summary_arr)[1] * dim(summary_arr)[2],
    ncol = dim(summary_arr)[3]
  )

  out_file <- file.path(out_dir, "posterior_R0_merged_chains")

  pdf(paste0(out_file, ".pdf"), width = 12, height = 10)
  plot(seq(1, T_total, length.out = T_total), apply(pooled_summary, 2, median), typ = "l",
      ylim = c(0, 10),
      xlab = "Days",
      ylab = "Basic reproduction number (R0)",
      main = paste0("Pooled R0 rate across chains (gamma = ", gamma_fixed, " )"), cex = 2)
  lines(1:T_total, apply(pooled_summary, 2, quantile, 0.025), col = "red", lty = 2)
  lines(1:T_total, apply(pooled_summary, 2, quantile, 0.975), col = "red", lty = 2)
  lines(1:T_total, exp(time_basis_0 %*% beta_generator) / 0.1, col = "blue", lwd = 2)
  # Legend to explain colors/linetypes
  legend("topright",
        legend = c("Posterior median", "95% CI", "True generator"),
        col    = c("black", "red", "blue"),
        lty    = c(1, 2, 1),
        bty    = "n")
  dev.off()
}