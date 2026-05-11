logdet_prox <- function(mat, eta) {
  eig <- eigen(mat)
  D <- eig$values
  Q <- eig$vectors
  return((eta/2) * (Q %*% diag((D + sqrt(D^2 + 4/eta))) %*% t(Q)))
}

soft_threshold <- function(X, tau, diag_penalty = FALSE) {
  X_new <- sign(X) * pmax(abs(X) - tau, 0)
  if (!diag_penalty) diag(X_new) <- diag(X)
  return(X_new)
}

squared_norm <- function(mat) {
  if (length(dim(mat)) == 3) {
    n_tp <- dim(mat)[3]
    snorm <- sum(sapply(1:n_tp, function(t) norm(mat[,,t], type = "f")^2))
  }
  else {
    snorm <- norm(mat, type = "f")^2
  }
  return(snorm)
}

update_rho <- function(rho, rnorm, snorm, mu = 10, tau_inc = 2, tau_dec = 2) {
  if (rnorm > mu * snorm)
    return (tau_inc * rho)
  else if (snorm > mu * rnorm)
    return (rho / tau_dec)
  return (rho)
}

tvgl_R <- function(S, n, lambda = 0, beta = 0, penalty_type = "l1", rho = 1, tol = 1e-4, rtol = 1e-4,
                   max_iter = 1000, verbose = 0) {
  start_time <- Sys.time()
  ### Initialize variables
  # The number of time points
  n_tp <- dim(S)[3]
  # The number of variables
  p <- dim(S)[1]
  Theta <- array(0, dim = c(p, p, n_tp))
  Z0 <- Z0_old <- U0 <- array(0, dim = c(p, p, n_tp))
  Z1 <- Z1_old <- Z2 <- Z2_old <- U1 <- U2 <- array(0, dim = c(p, p, n_tp-1))
  divider <- c(2, rep(3, n_tp - 2), 2)
  for (iter in 1:max_iter) {
    A <- Z0 - U0
    A[,,1:(n_tp-1)] <- A[,,1:(n_tp-1)] + Z1 - U1
    A[,,2:n_tp] <- A[,,2:n_tp] + Z2 - U2
    eta <- n / (divider * rho)
    A <- sapply(1:n_tp, function(t) A[,,t] / divider[t], simplify = "array")
    A <- (A + sapply(1:n_tp, function(t) t(A[,,t]), simplify = "array")) / 2
    A <- sapply(1:n_tp, function(t) A[,,t] * (1 / eta[t]), simplify = "array")
    A <- A - S
    Theta <- sapply(1:n_tp, function(t) logdet_prox(A[,,t], eta[t]), simplify = "array")
    # Update Z0 (sparsness).
    Z0 <- sapply(1:n_tp, function(t) soft_threshold(Theta[,,t] + U0[,,t], lambda / rho, diag_penalty = FALSE), simplify = "array")
    # Update Z1 and Z2 (penalty over the time).
    A1 <- Theta[,,1:(n_tp-1)] + U1
    A2 <- Theta[,,2:n_tp] + U2
    if (penalty_type == "l1") {
      E <- sapply(1:(n_tp-1), function(t) soft_threshold(A2[,,t] - A1[,,t], 2 * beta / rho, diag_penalty = TRUE), simplify = "array")
    }
    else if (penalty_type == "laplacian") {
      E <- (A2 - A1) / (1 + 2 * (2 * beta / rho))
    }
    Z1 <- (A1 + A2 - E) / 2
    Z2 <- (A1 + A2 + E) / 2
    # Update U0, U1 and U2.
    U0 <- U0 + Theta - Z0
    U1 <- U1 + Theta[,,1:(n_tp-1)] - Z1
    U2 <- U2 + Theta[,,2:n_tp] - Z2
    # Calculate norms.
    rnorm <- sqrt(squared_norm(Theta - Z0) + squared_norm(Theta[,,1:(n_tp-1)] - Z1) + squared_norm(Theta[,,2:n_tp] - Z2))
    snorm <- rho * sqrt(squared_norm(Z0 - Z0_old) + squared_norm(Z1 - Z1_old) + squared_norm(Z2 - Z2_old))
    Z0_old <- Z0
    Z1_old <- Z1
    Z2_old <- Z2
    # Calculate residuals.
    e_pri <- sqrt(prod(dim(Theta)) + 2 * prod(dim(Z1))) * tol + rtol * max(sqrt(squared_norm(Z0) + squared_norm(Z1) + squared_norm(Z2)),
                                                                           sqrt(squared_norm(Theta) + squared_norm(Theta[,,1:(n_tp-1)]) + squared_norm(Theta[,,2:n_tp])))
    e_dual <- sqrt(prod(dim(Theta)) + 2 * prod(dim(Z1))) * tol + rtol * rho * sqrt(squared_norm(U0) + squared_norm(U1) + squared_norm(U2))
    # Print diagnostics.
    if (verbose == 1 & iter %% 10 == 0) {
      lap_time <- Sys.time()
      elap_time <- as.numeric(difftime(lap_time, start_time, unit = "s"))
      cat("Iteration: ", iter, ". Elapsed time: ", round(elap_time, 3), " s. rnorm = ", round(rnorm, 3), ", e_pri = ", round(e_pri, 3), ", snorm = ", round(snorm, 3), ", e_dual = ", round(e_dual, 3), ", rho = ", rho, "\n", sep = "")
    }
    else if (verbose > 1) {
      lap_time <- Sys.time()
      elap_time <- as.numeric(difftime(lap_time, start_time, unit = "s"))
      cat("Iteration: ", iter, ". Elapsed time: ", round(elap_time, 3), " s. rnorm = ", round(rnorm, 3), ", e_pri = ", round(e_pri, 3), ", snorm = ", round(snorm, 3), ", e_dual = ", round(e_dual, 3), ", rho = ", rho, "\n", sep = "")
    }
    # Check convergence.
    if (rnorm <= e_pri & snorm <= e_dual) {
      end_time <- Sys.time()
      total_time <- as.numeric(difftime(end_time, start_time, unit = "s"))
      if (verbose >= 0) cat("Convergence reached at iteration ", iter, ". Elapsed time ", round(total_time, 3), " s.\n", sep = "")
      break
    }
    # Update rho.
    rho_new <- update_rho(rho, rnorm, snorm)
    scale <- rho / rho_new
    rho <- rho_new
    U0 <- U0 * scale
    U1 <- U1 * scale
    U2 <- U2 * scale
  }
  return_list <- list(Omega_ests = Z0, iters = iter, total_time = total_time)
  return(return_list)
}