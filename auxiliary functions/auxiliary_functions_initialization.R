# ──────────────────────────────────────────────────────────────────────────────
# Auxiliary functions
# ──────────────────────────────────────────────────────────────────────────────

# A simple jitter‐and‐select function
jitter_init <- function(map_pars, jitter_frac = 0.10, noise = "gaussian") {
  # pick the beta_free vector
  init_pars <- map_pars["beta_free"][[1]]
  # jitter each element relative to its own magnitude
  apply(init_pars, 1, function(x) {
    sd_vec <- pmax(abs(x) * jitter_frac, jitter_frac * 1e-1)
    if (length(x) == 1L) {
      if (noise == "gaussian") {
        return(rnorm(1, mean = x, sd = sd_vec))
      }
      # else, use uniform noise
      return(x + runif(1, -sd_vec, sd_vec))
      rnorm(1, mean = x, sd = sd_vec)
    } else {
      x + rnorm(length(x), mean = 0, sd = sd_vec)
    }
  })
}
jitter_init_2 <- function(x, jitter_frac = 0.10, noise = "gaussian") {
  # jitter each element relative to its own magnitude
  u_x <- log(x)
  sd_vec <- pmax(abs(u_x) * jitter_frac, jitter_frac * 1e-1)
  if (length(u_x) == 1L) {
    if (noise == "gaussian") {
      return(exp(rnorm(1, mean = u_x, sd = sd_vec)))
    }
    # else, use uniform noise
    return(u_x + runif(1, -sd_vec, sd_vec))
  } else {
    exp(u_x + rnorm(length(u_x), mean = 0, sd = sd_vec))
  }
}

generate_lhs_inits <- function(M, r, Sigma, theta_map) {
  # Args:
  #   M         : integer, total number of samples (must be S * r)
  #   r         : integer, reps per stratum
  #   Sigma     : p x p covariance matrix (positive definite)
  #   theta_map : length-p numeric vector (MAP in unconstrained space)
  #
  # Returns:
  #   M x p matrix of draws from N(theta_map, Sigma) via LHS with r reps per stratum
  
  p <- length(theta_map)
  if (M %% r != 0) {
      stop("M must be a multiple of r so that S = M/r is integer.")
  }
  S <- M / r
  
  # 1. Cholesky factor: Sigma = L %*% t(L)
  R <- chol(Sigma)     
  L <- t(R)
  
  # 2. Build the LHS uniform [0,1] matrix U (M x p)
  U <- matrix(NA_real_, nrow = M, ncol = p)
  for (i in seq_len(p)) {
      # stratify [0,1] into S intervals:
      cuts <- seq(0, 1, length.out = S + 1)
      # within each stratum, draw r uniforms
      u_strata <- unlist(lapply(seq_len(S), function(s) {
      runif(r, min = cuts[s], max = cuts[s+1])
      }))
      # randomly permute the M = S*r draws
      U[, i] <- sample(u_strata, size = M, replace = FALSE)
  }
  
  # 3. Transform to standard normal quantiles
  Z <- qnorm(U)          
  # 4. Map back to theta-space: theta = theta_map + L %*% z
  Theta <- sweep(Z %*% t(L), 2, theta_map, FUN = "+")
  
  colnames(Theta) <- names(theta_map)
  return(Theta)
}

rdirichlet0 <- function(n, alpha) {
  k <- length(alpha)
  t(sapply(1:n, function(i) {
    g <- rgamma(k, shape = alpha, rate = 1)
    g / sum(g)
  }))
}

# A helper that forks one Stan run, then waits up to `timeout` seconds:
run_with_timeout <- function(model, data, init, timeout = 5 * 60, returning = "value") {
  # Spawn the child
  job <- parallel::mcparallel({
    rstan::optimizing(
        model,
        data      = data,
        init      = init,
        verbose   = TRUE,
        tol_rel_obj = 0.5e-4,
        tol_obj     = 0.1,
        as_vector = FALSE
      )
  })
  
  # Wait for it up to `timeout` seconds
  res_list <- parallel::mccollect(job, wait = FALSE, timeout = timeout)
  
  if (is.null(res_list[[1]])) {
    # Timed out
    message("Job exceeded ", timeout, " seconds; killing fork.")
    tools::pskill(job$pid, tools::SIGKILL)
    return(NULL)
  }
  
  # Otherwise, return the optimization result
  if(returning == "value") {
    return(res_list[[1]]$value)
  } else{
    return(res_list[[1]])
  }
}

softmax <- function(x){
  exp(x) / sum(exp(x))
}

#' Hierarchical sampler for a TV–ball inside the simplex
#'
#' Generates `N_samples` probability vectors q that lie in the
#' total-variation ball of radius `dtv` centred at `p_ref`, i.e.
#'   ‖q − p_ref‖₁ ≤ 2 · dtv
#' and satisfy the simplex constraints (q ≥ 0, ∑q = 1).
#'
#' **Algorithm (greedy / sequential budget consumption)**
#'   1. Start with the full L¹-budget `rem ← 2 · dtv`.
#'   2. For the first *K − 1* coordinates draw
#'      \eqn{q_i \sim \mathrm{Unif}\!\bigl[\max(p_i − rem, 0),
#'      \min(p_i + rem, 1)\bigr]}.
#'      After each draw reduce the remaining budget by
#'      `abs(p_ref[i] − q[i])`.
#'   3. Set the last coordinate so the components sum to 1.
#'
#' The sampler is fast and always feasible, but **not uniform** in the ball.
#'
#' @param p_ref     Numeric vector of reference probabilities (must sum to 1).
#' @param dtv       Non-negative total-variation radius.
#' @param N_samples Integer ≥ 1; number of draws to return.
#'
#' @return A matrix with `N_samples` rows (or a plain vector if
#'         `N_samples == 1`).
#' @examples
#' set.seed(42)
#' hierarchical_DTV_ball_sampler(rep(1/3, 3), dtv = 0.15, N_samples = 5)
# Hierarchical sampler for q with ||q - p_ref||_1 <= 2*dtv and sum(q)=1
# Not uniform; simple feasibility sampler, collects up to N_Samples draws.

hierarchical_DTV_ball_sampler <- function(p_ref, dtv, N_Samples = 1L, max_trials = 100000L) {
  nq <- length(p_ref)
  if (nq < 2L) stop("p_ref must have at least two components")
  if (dtv < 0)  stop("dtv must be non-negative")
  if (N_Samples < 1L) stop("N_Samples must be >= 1")

  out <- matrix(NA_real_, nrow = N_Samples, ncol = nq)
  got <- 0L

  for (trial in seq_len(max_trials)) {
    q   <- numeric(nq)
    rem <- 2 * dtv  # TV radius converted to L1 budget

    valid <- TRUE
    # sample first nq-1 entries via deltas
    for (i in seq_len(nq - 1L)) {
      delta <- runif(1L, -rem, rem)
      qi    <- p_ref[i] + delta
      if (qi < 0 || qi > 1) { valid <- FALSE; break }  # invalid → retry
      q[i]  <- qi
      rem   <- rem - abs(delta)  # spend budget
    }

    if (valid) {
      # last coordinate: close the simplex
      q[nq] <- 1 - sum(q[-nq])

      # accept if inside [0,1] and last deviation fits remaining budget
      if (q[nq] >= 0 && q[nq] <= 1 && abs(q[nq] - p_ref[nq]) <= rem + 1e-12) {
        got <- got + 1L
        out[got, ] <- q
        if (got == N_Samples) break  # collected enough → stop early
      }
    }
  }

  if (got < N_Samples) {
    stop(sprintf("Sampling failed: obtained %d/%d samples within max_trials=%d",
                 got, N_Samples, max_trials))
  }

  if (N_Samples == 1L) drop(out) else out
}

build_jittered_init <- function(i, u_map, jittered_props, alpha_fixed_flag = FALSE, jitter_percentage,
                                gamma_fixed_flag = TRUE,
                                initialize_epidem = TRUE,
                                alpha_inits, gamma_inits,
                                jitter_init_epidem = TRUE,
                                jitter_init_props = TRUE
                                ) {
  li <- list(
    beta_free = jitter_init(u_map,
                            jitter_frac = jitter_percentage,
                            noise = "gaussian"))

  # θ
  if (jitter_init_props) {
    li$theta_free <- as.numeric(jittered_props[i, ])
  }

  # epidemic params
  if (initialize_epidem) {
    li$alpha_free <- if (alpha_fixed_flag) alpha_inits[1] else
      if (jitter_init_epidem) jitter_init_2(u_map$alpha_free,
                                            jitter_frac = jitter_percentage,
                                            noise = "gaussian") else u_map$alpha_free
    li$gamma_free <- if (gamma_fixed_flag) gamma_inits[1] else
      if (jitter_init_epidem) jitter_init_2(u_map$gamma_free,
                                            jitter_frac = jitter_percentage,
                                            noise = "gaussian") else u_map$gamma_free
    li$phi_inv_free <- if (jitter_init_epidem)
      jitter_init_2(u_map$phi_inv_free,
                    jitter_frac = jitter_percentage,
                    noise = "gaussian") else u_map$phi_inv_free
  }
  li$tau <- u_map$tau
  li
}

safe_value <- function(res) {

  if (is.list(res)) {
    v <- res[["value"]]              # [[ ]] won’t partial‑match
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) return(v)
  }

  -Inf
}