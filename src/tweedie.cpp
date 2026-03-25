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

arma::vec get_log_W(const arma::vec& alpha, int j, const arma::vec& constant_log_W) {
  return j * (constant_log_W - (1 + alpha) * log(j)) - log(2 * M_PI) -
         0.5 * arma::log(alpha) - log(j);
}

// Compute the logarithm of the whole sum
arma::vec log_A(const arma::vec& y, const arma::vec& phi, const arma::vec& rho) {
  arma::vec alpha = (2 - rho) / (rho - 1);
  arma::vec log_z = alpha % arma::log(y) - alpha % arma::log(rho - 1) - arma::log(2 - rho) -
                    (1 + alpha) % arma::log(phi);

  arma::vec j_max = arma::pow(y, 2 - rho) / (phi % (2 - rho));
  arma::vec log_W_max =
      j_max % (1 + alpha) - log(2 * M_PI) - 0.5 * arma::log(alpha) - arma::log(j_max);

  arma::vec constant_log_W = log_z + (1 + alpha) - alpha % arma::log(alpha);

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
  arma::mat log_W = j * log_z.t();
  log_W.each_col() -= arma::lgamma(j + 1);
  log_W -= arma::lgamma(j * alpha.t());

  arma::rowvec max_log_W = arma::max(log_W, 0);
  arma::mat exp_stabilized = arma::exp(log_W.each_row() - max_log_W);
  arma::rowvec log_sum_W = max_log_W + arma::log(arma::sum(exp_stabilized, 0));

  return log_sum_W.t() - arma::log(y);
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

  arma::uvec neg_idx = arma::find(x < 0);
  arma::uvec zero_idx = arma::find(x == 0);
  arma::uvec pos_idx = arma::find(x > 0);

  if (!neg_idx.is_empty()) {
    log_p(neg_idx).fill(-arma::datum::inf);
  }

  if (!zero_idx.is_empty()) {
    arma::vec mean_zero = mean(zero_idx);
    arma::vec power_zero = power(zero_idx);
    arma::vec dispersion_zero = dispersion(zero_idx);

    log_p(zero_idx) =
        -(arma::pow(mean_zero, 2 - power_zero)) / (dispersion_zero % (2 - power_zero));
  }

  if (!pos_idx.is_empty()) {
    arma::vec mean_pos = mean(pos_idx);
    arma::vec power_pos = power(pos_idx);
    arma::vec dispersion_pos = dispersion(pos_idx);

    log_p(pos_idx) = log_A(x(pos_idx), dispersion_pos, power_pos) +
                      ((x(pos_idx) % (arma::pow(mean_pos, 1 - power_pos) / (1 - power_pos))) -
                       (arma::pow(mean_pos, 2 - power_pos)) / (2 - power_pos)) /
                          dispersion_pos;
  }

  arma::vec result = log ? log_p : arma::exp(log_p);
  return result;
}




// Fully evaluate the cumulative density function
arma::vec dpois(const arma::vec& k, const arma::vec& lambda, bool log_p = false) {
  arma::vec out(k.n_elem, arma::fill::none);
  for (arma::uword i = 0; i < k.n_elem; ++i) {
    out[i] = R::dpois(k[i], lambda[i], log_p);
  }
  return out;
}

arma::vec pgamma(const arma::vec& x, const arma::vec& shape, const arma::vec& beta) {
  arma::vec out(x.n_elem, arma::fill::none);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    out[i] = R::pgamma(x[i], shape[i], 1.0 / beta[i], true, false);
  }
  return out;
}

arma::vec compound_Poisson_Gamma(const arma::vec& x,
                                 const arma::vec& lambda,
                                 const arma::vec& alpha,
                                 const arma::vec& beta) {
  arma::uword n = x.n_elem;
  arma::vec cdf = arma::exp(-lambda); 

  arma::vec l_mode = arma::floor(lambda);
  l_mode.elem(arma::find(l_mode < 1)).fill(1);

  arma::vec log_p_mode = dpois(l_mode, lambda, true);
  cdf += arma::exp(log_p_mode) % pgamma(x, alpha % l_mode, beta);

  arma::vec l_low = l_mode - 1;
  arma::uvec idx_low = arma::find(l_low > 0);
  while (!idx_low.is_empty()) {
    arma::vec l_low_active = l_low.elem(idx_low);
    arma::vec log_p_low = dpois(l_low_active, lambda.elem(idx_low), true);
    arma::vec keep_metric = log_p_mode.elem(idx_low) - log_p_low;
    arma::uvec keep = arma::find(keep_metric < 37.0);
    if (keep.is_empty()) {
      break;
    }

    arma::uvec active = idx_low.elem(keep);
    arma::vec l_active = l_low.elem(active);
    cdf.elem(active) += arma::exp(log_p_low.elem(keep)) %
      pgamma(x.elem(active), alpha.elem(active) % l_active, beta.elem(active));

    l_low.elem(active) -= 1;
    idx_low = arma::find(l_low > 0);
  }

  arma::vec l_high = l_mode + 1;
  arma::uvec idx_high = arma::regspace<arma::uvec>(0, n - 1);
  while (!idx_high.is_empty()) {
    arma::vec l_high_active = l_high.elem(idx_high);
    arma::vec log_p_high = dpois(l_high_active, lambda.elem(idx_high), true);
    arma::vec keep_metric = log_p_mode.elem(idx_high) - log_p_high;
    arma::uvec keep = arma::find(keep_metric < 37.0);
    if (keep.is_empty()) {
      break;
    }

    arma::uvec active = idx_high.elem(keep);
    arma::vec l_active = l_high.elem(active);
    cdf.elem(active) += arma::exp(log_p_high.elem(keep)) %
      pgamma(x.elem(active), alpha.elem(active) % l_active, beta.elem(active));

    l_high.elem(active) += 1;
    idx_high = active;
  }

  return arma::clamp(cdf, 0.0, 1.0);
}

// [[Rcpp::export]]
arma::vec tweedieCDF(arma::vec x, arma::vec mean, arma::vec dispersion,
                    arma::vec power) {
  int l = std::max({static_cast<int>(x.n_elem), static_cast<int>(mean.n_elem),
                    static_cast<int>(dispersion.n_elem), static_cast<int>(power.n_elem)});
  x = recycle_to_length(x, l, "x");
  mean = recycle_to_length(mean, l, "mean");
  dispersion = recycle_to_length(dispersion, l, "dispersion");
  power = recycle_to_length(power, l, "power");

  arma::vec lambda = arma::pow(mean, 2 - power) / (dispersion % (2 - power));
  arma::vec alpha = (2 - power) / (power - 1);
  arma::vec beta = 1 / (dispersion % (power - 1) % arma::pow(mean, power - 1));
  
  arma::vec cdf(l, arma::fill::zeros);

  arma::uvec neg_idx = arma::find(x < 0);
  arma::uvec zero_idx = arma::find(x == 0);
  arma::uvec pos_idx = arma::find(x > 0);

  if (!neg_idx.is_empty()) {
    cdf(neg_idx).fill(0);
  }

  if (!zero_idx.is_empty()) {
    cdf(zero_idx) = arma::exp(-lambda(zero_idx));
  }

  if (!pos_idx.is_empty()) {
    cdf(pos_idx) = compound_Poisson_Gamma(x(pos_idx), lambda(pos_idx), 
                                          alpha(pos_idx), beta(pos_idx));
  }

  return cdf;
}