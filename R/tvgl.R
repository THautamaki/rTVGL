#' Time-Varying Graphical LASSO for R
#' 
#' This function performs the ADMM algorithm for time-varying graphical LASSO (TVGL or TVGLASSO)
#' introduced by Hallac et al. (2017) and adapted from the Python package \code{regain} 
#' (Tomasi et al., 2018).
#'
#' @param data A list or an array containing the data, where each list item or array "slice" is
#'   \code{n} by \code{p} matrix, where \code{n} = sample size and \code{p} = number of variables.
#'   If an array provided, time must be the last dimension of the array.
#' @param lambda Numeric, a regularization parameter that controls the sparsity of the networks.
#' @param beta Numeric, a regularization parameter that controls the how strongly adjacent time
#'   points are pulled together.
#' @param penalty_type What kind of penalization is used over the time. Possible choices are
#'   \code{"l1"} and \code{"laplacian"}. The default value is \code{"l1"}.
#' @param use_correlation If \code{TRUE}, the sample correlation is used instead of sample covariance
#'   and number of observations is set to be 1. The default value is \code{TRUE}.
#' @param rho An Augmented Lagrangian penalty parameter in the ADMM algorithm.
#' @param tol An absolute tolerance for the convergence. The default value is 1e-4.
#' @param rtol A relative tolerance for the convergence. The default value is 1e-4.
#' @param max_iterations A maximum numbers of the ADMM algorithm iterations. The default value is 500.
#' @param use_Cpp Boolean, which is used to select the C++ or R version of the algorithm. The default
#'   value is \code{TRUE} (C++ version).
#' @param verbose Numeric with possible values of 0, 1 and 2. 0 (default) prints only information of
#'   final iteration. 1 prints information every 10th iteration, and 2 prints information every
#'   iteration.
#'
#' @return
#' A list, which contains following objects:
#' \item{Omega_ests}{
#'    A \code{p} by \code{p} by \code{n_timepoints} array containing the precision matrix estimates.
#' }
#' \item{Theta_ests}{
#'    A \code{p} by \code{p} by \code{n_timepoints} array containing adjacency matrices of the
#'    network estimates.
#' }
#' \item{iters}{
#'    The number of ADMM algorithm iterations.
#' }
#' \item{total_time}{
#'    The total run time of the algorithm in seconds.
#' }
#' @export
#' @references Hallac, D., Park, Y., Boyd, S., & Leskovec, J. (2017). Network inference via the 
#'    time-varying graphical lasso. In \emph{Proceedings of the 23rd ACM SIGKDD international
#'    conference on knowledge discovery and data mining} (pp. 205-213).
#' 
#' Tomasi, F., Tozzo, V., Salzo, S., & Verri, A. (2018). Latent variable time-varying network
#'    inference. In \emph{Proceedings of the 24th ACM SIGKDD International Conference on Knowledge
#'    Discovery & Data Mining} (pp. 2338-2346).
#'
#' @examples
#' sim <- generate_timeseries_network_data(n = 50, p = 100, n_timepoints = 20,
#'                                         change_points = c(5, 10, 15), n_add_del = 10)
#' results <- tvgl(sim$datasets[[1]], lambda = 0.25, beta = 0.25, penalty_type = "l1")
tvgl <- function(data, lambda = 0, beta = 0, penalty_type = "l1", use_correlation = TRUE, rho = 1,
                 tol = 1e-4, rtol = 1e-4, max_iterations = 1000, use_Cpp = TRUE, verbose = 0) {
  # If dataset is in array format, convert to the list.
  if (is.array(data)) {
    n_tp <- dim(data)[3]
    data <- sapply(1:n_tp, function(t) data[,,t])
  }
  n <- sapply(data, nrow)
  # The number of time points.
  n_tp <- length(data)
  # The number of variables
  p <- ncol(data[[1]])
  # Calculate sample covariance (biased so that it can be calculated even with one observation).
  S <- sapply(1:n_tp, function (t) t(data[[t]]) %*% data[[t]] / n[t], simplify = "array")
  # Convert to the sample correlation, is use_correlation is TRUE and set sample size to 1.
  if (use_correlation) {
    S <- sapply(1:n_tp, function(t) cov2cor(S[,,t]), simplify = "array")
    n <- rep(1, n_tp)
  }
  # Run TVGLASSO algorithm.
  if (use_Cpp) {
    results <- tvgl_Cpp(S, n, lambda = lambda, beta = beta, penalty_type = penalty_type, rho = rho,
                        tol = tol, rtol = rtol, max_iter = max_iterations, verbose = verbose)
  }
  else {
    results <- tvgl_R(S, n, lambda = lambda, beta = beta, penalty_type = penalty_type, rho = rho,
                      tol = tol, rtol = rtol, max_iter = max_iterations, verbose = verbose)
  }
  # Calculate adjacency matrices.
  Theta_ests <- array(dim = c(p, p, n_tp))
  for (t in 1:n_tp) {
    for (t in 1:n_tp) {
      theta_est <- results$Omega_ests[,,t]
      theta_est[theta_est != 0] <- 1
      diag(theta_est) <- 0
      Theta_ests[,,t] <- theta_est
    }
  }
  results$Theta_ests <- Theta_ests
  return(results)
}