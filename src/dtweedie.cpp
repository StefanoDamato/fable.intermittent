#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

arma::vec recycle_to_length(const arma::vec& x, int l, const std::string& name) {
  if (x.n_elem == static_cast<arma::uword>(l)) {
    return x;
  }

  if (x.n_elem == 1) {
    arma::vec out(l);
    out.fill(x[0]);
    return out;
  }

  stop("`%s` must have length 1 or length %d.", name.c_str(), l);
}

// A function to compute an item of the sum
double get_log_W(double alpha, int j, double constant_log_W) {
  return j * (constant_log_W - (1 + alpha) * log(j)) - log(2 * M_PI) -
         0.5 * log(alpha) - log(j);
}

// Compute the logarithm of the whole sum
double log_A(double y, double phi, double rho) {
  double alpha = (2 - rho) / (rho - 1);
  double log_z = alpha * std::log(y) - alpha * log(rho - 1) - log(2 - rho) -
                 (1 + alpha) * log(phi);

  double j_max = std::pow(y, 2 - rho) / (phi * (2 - rho));
  double log_W_max =
      j_max * (1 + alpha) - log(2 * M_PI) - 0.5 * log(alpha) - std::log(j_max);

  double constant_log_W = log_z + (1 + alpha) - alpha * log(alpha);

  int j_U = std::max(1., ceil(j_max));
  double log_W_U = get_log_W(alpha, j_U, constant_log_W);
  while (log_W_max - log_W_U < 37) {
    j_U = j_U + 1;
    log_W_U = get_log_W(alpha, j_U, constant_log_W);
  }

  int j_L = std::max(1., floor(j_max));
  double log_W_L = get_log_W(alpha, j_L, constant_log_W);
  while ((log_W_max - log_W_L < 37) && (j_L > 1)) {
    j_L = j_L - 1;
    log_W_L = get_log_W(alpha, j_L, constant_log_W);
  }

  arma::vec j = arma::linspace(j_L, j_U, j_U - j_L + 1);
  arma::vec log_W = j * log_z - arma::lgamma(j + 1) - arma::lgamma(alpha * j);
  double max_log_W = arma::max(log_W);
  double diff_sum_W = arma::sum(arma::exp(log_W - max_log_W));
  double log_sum_W = max_log_W + std::log(diff_sum_W);

  return log_sum_W - std::log(y);
}

double log_p_nonzero_arma(double y, double mu, double phi,
                          double rho) {
  return log_A(y, phi, rho) + ((y * (std::pow(mu, 1 - rho) / (1 - rho))) -
                               (std::pow(mu, 2 - rho)) / (2 - rho)) /
                                  phi;
}

// Fully evaluate the density for zero and non-zero values
// [[Rcpp::export]]
arma::vec tweedieDensity(arma::vec x, arma::vec mean, arma::vec dispersion,
                         arma::vec power, bool log) {
  int l = std::max({static_cast<int>(x.n_elem), static_cast<int>(mean.n_elem),
                    static_cast<int>(dispersion.n_elem), static_cast<int>(power.n_elem)});
  x = recycle_to_length(x, l, "x");
  mean = recycle_to_length(mean, l, "mean");
  dispersion = recycle_to_length(dispersion, l, "dispersion");
  power = recycle_to_length(power, l, "power");

  arma::vec log_p(l, arma::fill::none);

  arma::uvec zero_idx = arma::find(x == 0);
  arma::uvec nonzero_idx = arma::find(x > 0);

  if (!zero_idx.is_empty()) {
    arma::vec mean_zero = mean(zero_idx);
    arma::vec power_zero = power(zero_idx);
    arma::vec dispersion_zero = dispersion(zero_idx);

    log_p(zero_idx) =
        -(arma::pow(mean_zero, 2 - power_zero)) / (dispersion_zero % (2 - power_zero));
  }

  if (!nonzero_idx.is_empty()) {
    for (arma::uword i = 0; i < nonzero_idx.n_elem; ++i) {
      arma::uword idx = nonzero_idx[i];
      log_p[idx] = log_p_nonzero_arma(x[idx], mean[idx], dispersion[idx], power[idx]);
    }
  }

  arma::vec result = log ? log_p : arma::exp(log_p);
  return result;
}


// // Fully evaluate the cumulative density function
// // [[Rcpp::export]]
// arma::vec tweedieCDF(arma::vec x, arma::vec mean, arma::vec dispersion,
//                     arma::vec power) {
//   int l = std::max({static_cast<int>(x.n_elem), static_cast<int>(mean.n_elem),
//                     static_cast<int>(dispersion.n_elem), static_cast<int>(power.n_elem)});
//   x = recycle_to_length(x, l, "x");
//   mean = recycle_to_length(mean, l, "mean");
//   dispersion = recycle_to_length(dispersion, l, "dispersion");
//   power = recycle_to_length(power, l, "power");
  
//   arma::vec cdf(l, arma::fill::zeros);
  
//   arma::vec lambda = arma::pow(mean, 2 - power) / (dispersion % (2 - power));
//   arma::vec alpha = (2 - power) / (power - 1);
//   arma::vec beta = 1 / (dispersion % (power - 1) % arma::pow(mean, power - 1));

//   arma::ivec center = arma::conv_to<arma::ivec>::from(arma::floor(lambda));
//   center.elem(arma::find(center < 1)).fill(1);
//   arma::vec log_p_center(l, arma::fill::none);
//   double log_ratio_tol = -37.0;

//   for (int i = 0; i < l; ++i) {
//     log_p_center[i] = R::dpois(static_cast<double>(center[i]), lambda[i], true);
//   }

//   int l_low = arma::min(center);
//   while (l_low > 0) {
//     arma::vec log_p_next(l, arma::fill::none);
//     for (int i = 0; i < l; ++i) {
//       log_p_next[i] =
//           R::dpois(static_cast<double>(l_low - 1), lambda[i], true);
//     }
//     if (arma::all(log_p_next - log_p_center < log_ratio_tol)) {
//       break;
//     }
//     --l_low;
//   }

//   int l_high = arma::max(center);
//   while (true) {
//     arma::vec log_p_next(l, arma::fill::none);
//     for (int i = 0; i < l; ++i) {
//       log_p_next[i] =
//           R::dpois(static_cast<double>(l_high + 1), lambda[i], true);
//     }
//     if (arma::all(log_p_next - log_p_center < log_ratio_tol)) {
//       break;
//     }
//     ++l_high;
//   }

//   int l_low_eff = std::max(1, l_low);

//   for (int i = 0; i < l; ++i) {
//     if (x[i] < 0) {
//       cdf[i] = 0.0;
//       continue;
//     }

//     double cdf_i = R::dpois(0.0, lambda[i], false);
//     for (int l_count = l_low_eff; l_count <= l_high; ++l_count) {
//       double p_pois = R::dpois(static_cast<double>(l_count), lambda[i], false);
//       double p_gamma = R::pgamma(x[i], l_count * alpha[i], 1.0 / beta[i], true, false);
//       cdf_i += p_pois * p_gamma;
//     }

//     cdf[i] = std::min(1.0, std::max(0.0, cdf_i));
//   }

//   return cdf;
// }