#' Time series network data generator
#' 
#' This function generates time series network data with a random network structure
#'
#' @param n The number of samples per time point. A scalar or a \code{n_timepoints}-long vector.
#' @param p The number of variables.
#' @param n_timepoints The number of time points.
#' @param n_edges The number of edges in the network. The default value is 0 meaning that number
#'   of edges is set to \code{p}.
#' @param n_datasets The number of generated datasets, if multiple replicates needed. The default
#'   value is 1.
#' @param change_points A vector containing time points where changes happen.
#' @param n_add_del The number of changing edges in each change point.
#' @param rho Initial value for the diagonal of the precision matrix.
#' @param w1 Lower bound for the weights used to determine value of the precision matrix.
#' @param w2 Upper bound for the weights used to determine value of the precision matrix.
#' @param smooth_precision If \code{TRUE}, generated precision matrices change smoothly over time.
#'   If \code{FALSE}, changes happen immediately in the change points. The default value is \code{TRUE}.
#' @param seed_net Seed for network structure generation.
#' @param seed_data Seed for datasets generation.
#'
#' @return
#' A list, which contains following objects:
#' \item{Sigmas}{
#'    A \code{p} by \code{p} by \code{n_timepoints} array containing simulated covariance matrices.
#' }
#' \item{Omegas}{
#'    A \code{p} by \code{p} by \code{n_timepoints} array containing simulated precision matrices.
#' }
#' \item{Thetas}{
#'    A \code{p} by \code{p} by \code{n_timepoints} array containing simulated adjacency matrices.
#' }
#' \item{datasets}{
#'    An \code{n_datasets}-long list, each containing \code{n_timepoints}-long list, where each
#'    list item is \code{n} by \code{p} matrix of the simulated data per time point.
#' }
#' @export
#'
#' @examples
#' sim <- generate_timeseries_network_data(n = 50, p = 100, n_timepoints = 20,
#'                                         change_points = c(5, 10, 15), n_add_del = 10)
generate_timeseries_network_data <- function(n = 100, p = 100, n_timepoints = 10, n_edges = 0,
                                             n_datasets = 1, change_points = c(5), n_add_del = 5,
                                             rho = 0.25, w1 = 0.1, w2 = 0.3, smooth_precision = TRUE,
                                             seed_net = 12345, seed_data = 123456) {
  # Check most critical values and initialise some values.
  set.seed(seed_net)
  if (length(n) == 1) {
    n <- rep(n, n_timepoints)
  }
  else if (length(n) != n_timepoints) {
    stop("Length of the sample size vector must be either 1 or equal to the number of time points.")
  }
  if (any(1 == change_points) | length(change_points) > (n_timepoints - 1)) {
    stop("The first time point cannot be change point or too many change points.")
  }
  if (any(1 == diff(change_points)) & smooth_precision) {
    stop("If smooth_precision = TRUE, distance between change points must be greater than 1 time point.")
  }
  if (n_edges == 0) {
    n_edges <- p
  }
  cp <- change_points
  n_tp <- n_timepoints
  # Calcuate distance between change points.
  n_n <- diff(c(cp, n_tp)) + 1
  # Draw initial values for edges.
  a <- runif(n_edges, w1, w2)
  I <- diag(1, p)
  Theta <- rho * I
  
  # Generate network for the first time point (actually full network at this point).
  G <- igraph::graph_from_adjacency_matrix(matrix(1, p, p), diag = FALSE, mode = "undirected")
  
  # Convert network to edge list (matrix which first two columns have node numbers and third column
  # has edge weights).
  E <- igraph::as_edgelist(G)
  
  # Draw sample from the rows of the edge list, to create network which has n_edges connections.
  ind <- sample(1:nrow(E), n_edges)
  Es <- E[ind, ]
  Theta[Es] <- Theta[Es] - a
  Theta[Es[, 2:1]] <- Theta[Es[, 2:1]] - a
  b <- diag(diag(Theta), p)
  diag(Theta) <- diag(Theta) + abs(rowSums(Theta - b))
  
  G_Theta <- igraph::graph_from_adjacency_matrix(Theta, diag = FALSE, mode = "undirected",
                                                 weighted = TRUE)
  E <- igraph::as_data_frame(G_Theta)
  colnames(E) <- c("from", "to", "t1")
  
  # Expand edge list to contains all edges except self loops (so it has also zeros, but not diagonal
  # in adjacency matrix representation).
  E2 <- data.frame(from = rep(1:p, rep_len(p, p)), to = rep(1:p, p), t1 = rep(0, p*p))
  for (i in 1:nrow(E)) {
    inds <- E[i,1:2]
    E2[E2$from == inds$from & E2$to == inds$to, "t1"] <- E[i, "t1"]
  }
  E2 <- E2[which(E2$from != E2$to), 1:3]
  E <- E2
  
  # Expand edge list to contains all time points.
  E_temp <- matrix(0, nrow = nrow(E), ncol = n_tp - 1)
  colnames(E_temp) <- paste0("t", 2:n_tp)
  E <- cbind(E, E_temp)
  
  # Remove temp objects.
  rm(E_temp)
  rm(E2)
  
  # Set weights for the first time point.
  E[3:(cp[1] + 2)] = E[, 3]
  
  # If last time point is in the change points, for-loop goes over all change points.
  iterator <- 1:length(cp)
  # If last time point is not in the change points, add it and for-loop goes over change points
  # except that last time point.
  if (!any(n_tp == cp)) {
    cp <- c(cp, n_tp)
    iterator <- 1:(length(cp) - 1)
  }
  
  for (k in iterator) {
    # Check possible additions and deletions before any changes to avoid that same nodes to be
    # added and deleted at the same change points.
    add_ind <- which(E[, cp[k] + 2] == 0)
    del_ind <- which(E[, cp[k] + 2] != 0)
    # Draw sample for edges to be added in the network.
    add_edges <- sample(add_ind, n_add_del)
    # Create matrix which will have those edges.
    a_target_add <- matrix(0, n_n[k], n_add_del)
    # Draw weights for the edges.
    a_target_add[1, ] <- runif(n_add_del, min = w1, max = w2)
    # If smooth precision matrices needed, gradually increase the values. Otherwise, just replicate
    # values as many times as needed.
    if (smooth_precision) {
      a_target_add <- apply(a_target_add, 2, function(x) seq(0, x[1], length.out = n_n[k]))
    }
    else {
      a_target_add[1:n_n[k], ] <- t(replicate(n_n[k], a_target_add[1,]))
    }
    a_target_add <- t(a_target_add)
    # Set indices as row names.
    rownames(a_target_add) <- add_edges
    # Append values to the E matrix.
    if (cp[k] < n_tp) {
      E[add_edges, c(cp[k]:cp[k + 1] + 2)] <- -a_target_add
      E[-add_edges, c(cp[k]:cp[k + 1] + 2)] <- E[-add_edges, cp[k] + 2]
    }
    else {
      E[add_edges, cp[k] + 2] <- -a_target_add
      E[-add_edges, cp[k] + 2] <- E[-add_edges, cp[k] + 2]
    }
    # Delete edges from the network.
    del_edges <- sample(del_ind, n_add_del)
    a_target_del <- matrix(0, n_n[k], n_add_del)
    # If smooth precision matrices needed, gradually decrease the values.
    if (smooth_precision) {
      a_target_del[1, ] <- E[del_edges, cp[k] + 2]
      a_target_del <- apply(a_target_del, 2, function(x) seq(x[1], 0, length.out = n_n[k]))
    }
    a_target_del <- t(a_target_del)
    rownames(a_target_del) <- del_edges
    # Append values to the E matrix.
    if (cp[k] < n_tp) {
      E[del_edges, c(cp[k]:cp[k + 1] + 2)] <- a_target_del
      E[-c(del_edges, add_edges), c(cp[k]:cp[k + 1] + 2)] <- E[-c(del_edges, add_edges), cp[k] + 2]
    }
    else {
      E[del_edges, cp[k] + 2] <- a_target_del
      E[-c(del_edges, add_edges), cp[k] + 2] <- E[-c(del_edges, add_edges), cp[k] + 2]
    }
  }
  # Create covariance, precision and adjacency matrices.
  Omegas <- Sigmas <- Thetas <- array(NA, dim = c(p, p, n_tp))
  for (t in 1:n_tp) {
    E_temp <- as.matrix(E[, c(1, 2, t + 2)])
    G <- igraph::graph_from_edgelist(E_temp[, c(1, 2)])
    igraph::E(G)$weight <- E_temp[, 3]
    Omega <- igraph::as_adjacency_matrix(G, attr = "weight")
    Omega <- as.matrix(Omega)
    Omega <- (Omega + t(Omega))
    diag(Omega) <- rho
    b <- diag(diag(Omega), p)
    diag(Omega) <- diag(Omega) + abs(rowSums(Omega - b))
    Sigmas[,,t] <- cov2cor(solve(Omega))
    Omegas[,,t] <- solve(Sigmas[,,t])
    Theta <- Omega
    Theta[Theta != 0] <- 1
    diag(Theta) <- 0
    Thetas[,,t] <- Theta
  }
  # Generate datasets.
  datasets <- list()
  set.seed(seed_data)
  for (r in 1:n_datasets) {
    data <- list()
    for (t in 1:n_tp) {
      if (n[t] == 1) data[[t]] <- t(as.matrix(MASS::mvrnorm(n = n[t], mu = rep(0, p), Sigma = Sigmas[,,t])))
      else data[[t]] <- MASS::mvrnorm(n = n[t], mu = rep(0, p), Sigma = Sigmas[,,t])
    }
    datasets[[r]] <- data
  }
  return(list(Sigmas = Sigmas, Omegas = Omegas, Thetas = Thetas, datasets = datasets))
}
