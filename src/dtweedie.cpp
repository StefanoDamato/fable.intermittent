#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// A function to compute an item of the sum
arma::vec get_log_W(double alpha, int j, arma::vec constant_log_W) {
  return j * (constant_log_W - (1 + alpha) * log(j)) - log(2 * M_PI) -
         0.5 * log(alpha) - log(j);
}

// Compute the logarithm of the whole sum
arma::vec log_A(arma::vec y, double phi, double rho) {
  double alpha = (2 - rho) / (rho - 1);
  arma::vec log_z = alpha * arma::log(y) - alpha * log(rho - 1) - log(2 - rho) -
                    (1 + alpha) * log(phi);

  arma::vec j_max = arma::pow(y, 2 - rho) / (phi * (2 - rho));
  arma::vec log_W_max =
      j_max * (1 + alpha) - log(2 * M_PI) - 0.5 * log(alpha) - arma::log(j_max);

  arma::vec constant_log_W = log_z + (1 + alpha) - alpha * log(alpha);

  int j_U = std::max(1., ceil(arma::max(j_max)));
  arma::vec log_W_U = get_log_W(alpha, j_U, constant_log_W);
  while (any(log_W_max - log_W_U < 37)) {
    j_U = j_U + 1;
    log_W_U = get_log_W(alpha, j_U, constant_log_W);
  }

  int j_L = std::max(1., floor(arma::min(j_max)));
  arma::vec log_W_L = get_log_W(alpha, j_L, constant_log_W);
  while (any(log_W_max - log_W_L < 37) & (j_L > 1)) {
    j_L = j_L - 1;
    log_W_L = get_log_W(alpha, j_L, constant_log_W);
  }

  arma::vec j = arma::linspace(j_L, j_U, j_U - j_L + 1);
  arma::vec log_sum_W(log_z.n_elem);

  for (size_t i = 0; i < log_z.n_elem; i++) {
    arma::vec log_W =
        j * log_z[i] - arma::lgamma(j + 1) - arma::lgamma(alpha * j);
    double max_log_W = arma::max(log_W);
    double diff_sum_W = arma::sum(arma::exp(log_W - max_log_W));
    log_sum_W[i] = max_log_W + std::log(diff_sum_W);
  }

  return log_sum_W - arma::log(y);
}

arma::vec log_p_nonzero_arma(arma::vec y, arma::vec mu, double phi,
                             double rho) {
  return log_A(y, phi, rho) + ((y % (arma::pow(mu, 1 - rho) / (1 - rho))) -
                               (arma::pow(mu, 2 - rho)) / (2 - rho)) /
                                  phi;
}

// Fully evaluate the density for zero and non-zero values
// [[Rcpp::export]]
arma::vec tweedieDensity(arma::vec x, arma::vec mean, double dispersion,
                         double power, bool log) {
  int l = std::max(x.n_elem, mean.n_elem);
  x.resize(l);
  mean.resize(l);

  arma::vec log_p(l, arma::fill::none);

  arma::uvec zero_idx = arma::find(x == 0);
  arma::uvec nonzero_idx = arma::find(x > 0);

  if (!zero_idx.is_empty()) {
    log_p(zero_idx) =
        -(arma::pow(mean(zero_idx), 2 - power)) / (dispersion * (2 - power));
  }

  if (!nonzero_idx.is_empty()) {
    log_p(nonzero_idx) = log_p_nonzero_arma(x(nonzero_idx), mean(nonzero_idx),
                                            dispersion, power);
  }

  arma::vec result = log ? log_p : arma::exp(log_p);
  return result;
}
