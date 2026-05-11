tvgl <- function(data, lambda = 0, beta = 0, penalty_type = "l1", use_correlation = TRUE, rho = 1,
                 tol = 1e-4, rtol = 1e-4, max_iterations = 1000, use_Cpp = TRUE, verbose = 0) {
  # If dataset is in array format, convert to the list.
  if (is.array(data)) data <- lapply(1:n_tp, function(t) data[,,t])
  n <- sapply(data, nrow)
  # The number of time points.
  n_tp <- length(data)
  # The number of variables
  p <- ncol(data[[1]])
  # Calculate sample covariance (biased so that it can be calculated even with one observation).
  S <- sapply(1:n_tp, function (t) t(data[[t]]) %*% data[[t]] / n[t], simplify = "array")
  if (use_correlation) {
    S <- sapply(1:n_tp, function(t) cov2cor(S[,,t]), simplify = "array")
    n <- rep(1, n_tp)
  }
  if (use_Cpp) {
    results <- tvgl_Cpp(S, n, lambda = lambda, beta = beta, penalty_type = penalty_type, rho = rho,
                        tol = tol, rtol = rtol, max_iter = max_iterations, verbose = verbose)
  }
  else {
    results <- tvgl_R(S, n, lambda = lambda, beta = beta, penalty_type = penalty_type, rho = rho,
                      tol = tol, rtol = rtol, max_iter = max_iterations, verbose = verbose)
  }
  return(results)
}