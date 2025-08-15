# Load required packages
library(splines)

# Computes the penalty matrix K for the spline coefficients.
# Here, we use a difference penalty of order 'penalty_order':
#   K = t(D) %*% D, where D = diff(I_m, differences = penalty_order)
compute_K <- function(penalty_order, m) {
  if (m <= penalty_order) {
    stop("The number of basis functions 'm' must be greater than the penalty order.")
  }
  I_m <- diag(m)
  D <- diff(I_m, differences = penalty_order)
  K <- t(D) %*% D
  return(K)
}

# Computes the B-spline basis matrix evaluated at the given time points.
# Arguments:
#   times: a vector of time points (should be sorted)
#   bdegree: degree of the B-splines
#   m: total number of basis functions (degrees of freedom)

# Non necesarily equidistant knots
# compute_B_matrix <- function(times, bdegree, m) {
#   # Use the built-in bs() function from the splines package.
#   # Setting 'df = m' produces m basis functions.
#   B <- bs(times, degree = bdegree, df = m, intercept = TRUE)
#   return(as.matrix(B))
# }

# Equidistant knots
compute_B_matrix <- function(times, bdegree, m) {
  # Use the built-in bs() function from the splines package.
  # Setting 'df = m' produces m basis functions.
  # B <- bs(times, degree = bdegree, df = m, intercept = TRUE)
  # return(as.matrix(B))
  
  # New form
  # B-spline basis
  timesl <- min(times)
  timesr <- max(times)
  dtimes <- (timesr - timesl)/(m - bdegree)

  knots <- seq(timesl - bdegree*dtimes, timesr + bdegree*dtimes, by=dtimes)
  B <- splines::spline.des(knots, times, ord = bdegree + 1, outer.ok = TRUE)$design
  return(as.matrix(B))
}

# Convenience function that returns both the penalty matrix and the B-spline basis.
# Arguments:
#   times: a vector of time points (the grid on which the basis is evaluated)
#   bdegree: degree of the B-splines
#   m: total number of basis functions
#   penalty_order: order of the random walk (difference penalty)
compute_spline_objects <- function(times, bdegree, m, penalty_order) {
  K <- compute_K(penalty_order, m)
  B_matrix <- compute_B_matrix(times, bdegree, m)
  list(
    ts_bs = times,
    B_matrix = B_matrix,
    K = K
  )
}

# # Example usage:
# # Suppose we want to use cubic B-splines (degree 3) with m = 10 basis functions,
# # but we choose a random walk penalty of order 2.
# times_grid <- seq(0, 100, length.out = 50)  # e.g., 50 time points over [0, 100]
# bdegree <- 3       # B-spline degree
# m <- 10            # number of basis functions
# penalty_order <- 2 # order of the difference penalty

# spline_objs <- compute_spline_objects(times_grid, bdegree, m, penalty_order)
# print(spline_objs$K)              # roughness penalty matrix for Stan (K_bs)
# print(head(spline_objs$B_matrix))  # first few rows of the B-spline basis matrix (B)
