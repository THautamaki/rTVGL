#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace std::chrono;
using namespace arma;
using namespace Rcpp;

arma::mat logdet_prox(arma::mat matrix, double eta) {
  vec D;
  mat Q;
  eig_sym(D, Q, matrix);
  return (eta/2) * (Q * diagmat(D + sqrt(pow(D, 2) + 4/eta)) * trans(Q));
}

arma::mat soft_threshold(arma::mat matrix, double tau, bool diag_penalty = false) {
  mat matrix_new;
  int p = matrix.n_cols;
  matrix_new = sign(matrix) % max(abs(matrix) - tau, zeros(p, p));
  if (!diag_penalty) {
    matrix_new.diag() = matrix.diag();
  }
  return matrix_new;
}

double cube_squared_norm(arma::cube cube) {
  int n_tp = cube.n_slices;
  double sq_norm = 0;
  for (int t = 0; t < n_tp; t++) {
    sq_norm += pow(norm(cube.slice(t), "for"), 2);
  }
  return sq_norm;
}

double update_rho(double rho, double rnorm, double snorm, double mu = 10, double tau_inc = 2,
                  double tau_dec = 2) {
  if (rnorm > mu * snorm) {
    return tau_inc * rho;
  } else if (snorm > mu * rnorm) {
    return (rho / tau_dec);
  }
  return (rho);
}

// [[Rcpp::export]]
List tvgl_Cpp(const arma::cube& S, arma::vec n, double lambda = 0, double beta = 0,
              const std::string& penalty_type = "l1", double rho = 1,
              double tol = 1e-4, double rtol = 1e-4, int max_iter = 500, int verbose = 0) {
  auto start_time = high_resolution_clock::now();
  int n_tp = S.n_slices;
  int p = S.n_cols;
  cube Omega = zeros(p, p, n_tp);
  cube Z0 = zeros(p, p, n_tp);
  cube Z1 = zeros(p, p, n_tp-1);
  cube Z2 = zeros(p, p, n_tp-1);
  cube Z0_old = zeros(size(Z0));
  cube Z1_old = zeros(size(Z1));
  cube Z2_old = zeros(size(Z2));
  cube U0 = zeros(size(Z0));
  cube U1 = zeros(size(Z1));
  cube U2 = zeros(size(Z2));
  vec divider = ones(n_tp) * 3;
  divider(0) = 2;
  divider(n_tp - 1) = 2;
  int iters = 0;
  Rcout << fixed << setprecision(3);
  for (int i = 0; i < max_iter; i++) {
    iters += 1;
    cube A = Z0 - U0;
    A.slices(0, n_tp-2) += Z1 - U1;
    A.slices(1, n_tp-1) += Z2 - U2;
    vec eta = n / (divider * rho);
    for (unsigned int j = 0; j < A.n_slices; j++) {
      A.slice(j) /= divider(j);
    }
    for (unsigned int j = 0; j < A.n_slices; j++) {
      A.slice(j) += A.slice(j).t();
    }
    A /= 2;
    for (unsigned int j = 0; j < A.n_slices; j++) {
      A.slice(j) /= eta(j);
    }
    A -= S;
    Omega = A;
    for (unsigned int j = 0; j < Omega.n_slices; j++) {
      Omega.slice(j) = logdet_prox(Omega.slice(j), eta(j));
    }
    Z0 = Omega + U0;
    Z0 = sign(Z0) % max(abs(Z0) - lambda / rho, zeros(size(Z0)));
    for (unsigned int j = 0; j < Z0.n_slices; j++) {
      Z0.slice(j).diag() = Omega.slice(j).diag() + U0.slice(j).diag();
    }
    cube A1 = Omega.slices(0, n_tp-2) + U1;
    cube A2 = Omega.slices(1, n_tp-1) + U2;
    cube E = A2 - A1;
    if (penalty_type == "l1") {
      E =  sign(E) % max(abs(E) - 2 * beta / rho, zeros(size(E)));
    } else if (penalty_type == "laplacian") {
      E /= (1 + 2 * (2 * beta / rho));
    }
    Z1 = (A1 + A2 - E) / 2;
    Z2 = (A1 + A2 + E) / 2;
    U0 += Omega - Z0;
    U1 += Omega.slices(0, n_tp-2) - Z1;
    U2 += Omega.slices(1, n_tp-1) - Z2;
    double rnorm = sqrt(cube_squared_norm(Omega - Z0) + cube_squared_norm(Omega.slices(0, n_tp-2) - Z1) + cube_squared_norm(Omega.slices(1, n_tp-1) - Z2));
    double snorm = rho * sqrt(cube_squared_norm(Z0 - Z0_old) + cube_squared_norm(Z1 - Z1_old) + cube_squared_norm(Z2 - Z2_old));
    Z0_old = Z0;
    Z1_old = Z1;
    Z2_old = Z2;
    int size_Omega = Omega.n_rows * Omega.n_cols * Omega.n_slices;
    int size_Z1 = Z1.n_rows * Z1.n_cols * Z1.n_slices;
    double e_pri = sqrt(size_Omega + 2 * size_Z1) * tol + rtol * max(sqrt(cube_squared_norm(Z0) + cube_squared_norm(Z1) + cube_squared_norm(Z2)),
                        sqrt(cube_squared_norm(Omega) + cube_squared_norm(Omega.slices(0, n_tp-2)) + cube_squared_norm(Omega.slices(1, n_tp-1))));
    double e_dual = sqrt(size_Omega + 2 * size_Z1) * tol + rtol * rho * sqrt(cube_squared_norm(U0) + cube_squared_norm(U1) + cube_squared_norm(U2));
    if (verbose == 1 && iters % 10 == 0) {
      auto lap_time = high_resolution_clock::now();
      double elap_time = duration_cast<microseconds>(lap_time - start_time).count();
      Rcout << "Iteration: " << iters << ". Elapsed time: " << elap_time/1000000 << " s. rnorm = " <<
        rnorm << ", e_pri = " << e_pri << ", snorm = " << snorm << ", e_dual = " << e_dual << ", rho = " << rho << "\n";
    }
    else if (verbose == 2) {
      auto lap_time = high_resolution_clock::now();
      double elap_time = duration_cast<microseconds>(lap_time - start_time).count();
      Rcout << "Iteration: " << iters << ". Elapsed time: " << elap_time/1000000 << " s. rnorm = " <<
        rnorm << ", e_pri = " << e_pri << ", snorm = " << snorm << ", e_dual = " << e_dual << ", rho = " << rho << "\n";
    }
    if (rnorm <= e_pri && snorm <= e_dual) {
      auto end_time = high_resolution_clock::now();
      double elap_time = duration_cast<microseconds>(end_time - start_time).count();
      Rcout << "Convergence reached at iteration " << iters << ". Elapsed time " << elap_time/1000000 << " s.\n";
      break;
    }
    double rho_new = update_rho(rho, rnorm, snorm);
    double scale = rho / rho_new;
    rho = rho_new;
    U0 *= scale;
    U1 *= scale;
    U2 *= scale;
  }
  auto end_time = high_resolution_clock::now();
  double elap_time = duration_cast<microseconds>(end_time - start_time).count();
  List results = List::create(_["Omega_ests"] = Z0, _["iters"] = iters, _["tot_time"] = elap_time/1000000);
  return results;
}