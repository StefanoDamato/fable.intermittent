#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Helper function to reshape vectors
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

  // Set alpha and log(z) which will be used repeatedly
  arma::vec alpha = (2 - rho) / (rho - 1);
  arma::vec log_z = alpha % arma::log(y) - alpha % arma::log(rho - 1) - arma::log(2 - rho) -
                    (1 + alpha) % arma::log(phi);

  // Compute the location and the value of the peak term
  arma::vec j_max = arma::pow(y, 2 - rho) / (phi % (2 - rho));
  arma::vec log_W_max =
      j_max % (1 + alpha) - log(2 * M_PI) - 0.5 * arma::log(alpha) - arma::log(j_max);
  arma::vec constant_log_W = log_z + (1 + alpha) - alpha % arma::log(alpha);

  // Iterate to identify the lower bound
  int j_U = std::max(1., ceil(arma::max(j_max)));
  arma::vec log_W_U = get_log_W(alpha, j_U, constant_log_W);
  while (any(log_W_max - log_W_U < 37)) {
    j_U = j_U + 1;
    log_W_U = get_log_W(alpha, j_U, constant_log_W);
  }

  // Iterate to identify the upper bound
  int j_L = std::max(1., floor(arma::min(j_max)));
  arma::vec log_W_L = get_log_W(alpha, j_L, constant_log_W);
  while (any(log_W_max - log_W_L < 37) & (j_L > 1)) {
    j_L = j_L - 1;
    log_W_L = get_log_W(alpha, j_L, constant_log_W);
  }

  // Create the sum and evaluate the expression of log(W(j)) for each value
  arma::vec j = arma::linspace(j_L, j_U, j_U - j_L + 1);
  arma::mat log_W = j * log_z.t();
  log_W.each_col() -= arma::lgamma(j + 1);
  log_W -= arma::lgamma(j * alpha.t());

  // Divide by the maximum term and move to logarithm for stability
  arma::rowvec max_log_W = arma::max(log_W, 0);
  arma::mat exp_stabilized = arma::exp(log_W.each_row() - max_log_W);
  arma::rowvec log_sum_w = max_log_W + arma::log(arma::sum(exp_stabilized, 0));

  // Return the aggregated term
  return log_sum_w.t() - arma::log(y);
}


// Fully evaluate the density for zero and non-zero values
// [[Rcpp::export]]
arma::vec tweedieDensity(arma::vec x, arma::vec mean, arma::vec dispersion,
                         arma::vec power, bool log) {

  // Make all inputs of the same length
  int l = std::max({static_cast<int>(x.n_elem), static_cast<int>(mean.n_elem),
                    static_cast<int>(dispersion.n_elem), static_cast<int>(power.n_elem)});
  x = recycle_to_length(x, l, "x");
  mean = recycle_to_length(mean, l, "mean");
  dispersion = recycle_to_length(dispersion, l, "dispersion");
  power = recycle_to_length(power, l, "power");
  arma::vec log_p(l, arma::fill::none);

  // For negatives the density is 0
  arma::uvec neg_idx = arma::find(x < 0);
  if (!neg_idx.is_empty()) {
    log_p(neg_idx).fill(-arma::datum::inf);
  }

  // In 0 the density is exp(-lambda)
  arma::uvec zero_idx = arma::find(x == 0);
  if (!zero_idx.is_empty()) {
    arma::vec mean_zero = mean(zero_idx);
    arma::vec power_zero = power(zero_idx);
    arma::vec dispersion_zero = dispersion(zero_idx);
    log_p(zero_idx) =
        -(arma::pow(mean_zero, 2 - power_zero)) / (dispersion_zero % (2 - power_zero));
  }

  // On the positive axis evaluate the density and the sum in A
  arma::uvec pos_idx = arma::find(x > 0);
  if (!pos_idx.is_empty()) {
    arma::vec mean_pos = mean(pos_idx);
    arma::vec power_pos = power(pos_idx);
    arma::vec dispersion_pos = dispersion(pos_idx);
    log_p(pos_idx) = log_A(x(pos_idx), dispersion_pos, power_pos) +
                      ((x(pos_idx) % (arma::pow(mean_pos, 1 - power_pos) / (1 - power_pos))) -
                       (arma::pow(mean_pos, 2 - power_pos)) / (2 - power_pos)) /
                          dispersion_pos;
  }

  // Return the final vector with log-trandormation if needed
  arma::vec result = log ? log_p : arma::exp(log_p);
  return result;
}


// Vectorise the Poisson mass function
arma::vec dpois(const arma::vec& k, const arma::vec& lambda, bool log_p = false) {
  arma::vec out(k.n_elem, arma::fill::none);
  for (arma::uword i = 0; i < k.n_elem; ++i) {
    out[i] = R::dpois(k[i], lambda[i], log_p);
  }
  return out;
}


// Vectorise the Gamma density function
arma::vec pgamma(const arma::vec& x, const arma::vec& shape, const arma::vec& beta) {
  arma::vec out(x.n_elem, arma::fill::none);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    out[i] = R::pgamma(x[i], shape[i], 1.0 / beta[i], true, false);
  }
  return out;
}

// Evaluate the compound Poisson-Gamma cumulative density with truncation
arma::vec compound_Poisson_Gamma(const arma::vec& x,
                                 const arma::vec& lambda,
                                 const arma::vec& alpha,
                                 const arma::vec& beta) {

  // Identify the mode of the Poisson                                
  arma::uword n = x.n_elem;
  arma::vec l_mode = arma::floor(lambda);
  l_mode.elem(arma::find(l_mode < 1)).fill(1);
  arma::vec log_p_mode = dpois(l_mode, lambda, true);
  
  // Evaluate the CDF at the mode and add the mass in 0
  arma::vec cdf = arma::exp(-lambda); 
  cdf += arma::exp(log_p_mode) % pgamma(x, alpha % l_mode, beta);

  // Iterate on the values of the Poisson below the mode
  arma::vec l_low = l_mode - 1;
  arma::uvec idx_low = arma::find(l_low > 0);
  while (!idx_low.is_empty()) {

    // Check for which indices the summand are not negligible
    arma::vec l_low_active = l_low.elem(idx_low);
    arma::vec log_p_low = dpois(l_low_active, lambda.elem(idx_low), true);
    arma::vec keep_metric = log_p_mode.elem(idx_low) - log_p_low;
    arma::uvec keep = arma::find(keep_metric < 37.0);
    if (keep.is_empty()) {
      break;
    }

    // For active indices evaluate the term of the CDF and add it
    arma::uvec active = idx_low.elem(keep);
    arma::vec l_active = l_low.elem(active);
    cdf.elem(active) += arma::exp(log_p_low.elem(keep)) %
      pgamma(x.elem(active), alpha.elem(active) % l_active, beta.elem(active));
    l_low.elem(active) -= 1;
    idx_low = arma::find(l_low > 0);
  }

  // Iterate on the values of the Poisson below the mode
  arma::vec l_high = l_mode + 1;
  arma::uvec idx_high = arma::regspace<arma::uvec>(0, n - 1);
  while (!idx_high.is_empty()) {

    // Check for which indices the summand are not negligible
    arma::vec l_high_active = l_high.elem(idx_high);
    arma::vec log_p_high = dpois(l_high_active, lambda.elem(idx_high), true);
    arma::vec keep_metric = log_p_mode.elem(idx_high) - log_p_high;
    arma::uvec keep = arma::find(keep_metric < 37.0);
    if (keep.is_empty()) {
      break;
    }

    // For active indices evaluate the term of the CDF and add it
    arma::uvec active = idx_high.elem(keep);
    arma::vec l_active = l_high.elem(active);
    cdf.elem(active) += arma::exp(log_p_high.elem(keep)) %
      pgamma(x.elem(active), alpha.elem(active) % l_active, beta.elem(active));
    l_high.elem(active) += 1;
    idx_high = active;
  }

  // Return the cdf with numerical correction
  return arma::clamp(cdf, 0.0, 1.0);
}

// [[Rcpp::export]]
arma::vec tweedieCDF(arma::vec x, arma::vec mean, arma::vec dispersion,
                    arma::vec power) {

  // Make all vector of the same size
  int l = std::max({static_cast<int>(x.n_elem), static_cast<int>(mean.n_elem),
                    static_cast<int>(dispersion.n_elem), static_cast<int>(power.n_elem)});
  x = recycle_to_length(x, l, "x");
  mean = recycle_to_length(mean, l, "mean");
  dispersion = recycle_to_length(dispersion, l, "dispersion");
  power = recycle_to_length(power, l, "power");

  // Compute the alternative parametrization
  arma::vec lambda = arma::pow(mean, 2 - power) / (dispersion % (2 - power));
  arma::vec alpha = (2 - power) / (power - 1);
  arma::vec beta = 1 / (dispersion % (power - 1) % arma::pow(mean, power - 1));
  arma::vec cdf(l, arma::fill::zeros);

  // For negatives the CDF is 0
  arma::uvec neg_idx = arma::find(x < 0);
  if (!neg_idx.is_empty()) {
    cdf(neg_idx).fill(0);
  }

  // In 0 the CDF is exp(-lambda)
  arma::uvec zero_idx = arma::find(x == 0);
  if (!zero_idx.is_empty()) {
    cdf(zero_idx) = arma::exp(-lambda(zero_idx));
  }

  // For positive values use the compund Gamma-Poisson structure with truncation
  arma::uvec pos_idx = arma::find(x > 0);
  if (!pos_idx.is_empty()) {
    cdf(pos_idx) = compound_Poisson_Gamma(x(pos_idx), lambda(pos_idx), 
                                          alpha(pos_idx), beta(pos_idx));
  }

  // Return the value
  return cdf;
}


// Vectorized bisection root-finder for monotone CDF inversion
arma::vec solver_bisection(const arma::vec& q, const arma::vec& lambda, const arma::vec& alpha,
                          const arma::vec& beta, const arma::vec& x0) {

  // Specify some parameters
  int n = q.n_elem;
  double tol = 1e-8;
  int max_iter = 100;

  // Initialize the bounds in the same points
  arma::vec f_L = compound_Poisson_Gamma(x0, lambda, alpha, beta) - q;
  arma::vec f_U = f_L;
  arma::vec x_L = x0;
  arma::vec x_U = x0;

  // Expand the upper bound until f(x_U) >= 0 for all
  arma::uvec neg_idx = arma::find(f_U < 0);
  while (!neg_idx.is_empty()) {  
    x_L.elem(neg_idx) = x_U.elem(neg_idx);
    f_L.elem(neg_idx) = f_U.elem(neg_idx);
    x_U.elem(neg_idx) *= 2;
    f_U.elem(neg_idx) = compound_Poisson_Gamma(x_U.elem(neg_idx), lambda.elem(neg_idx), alpha.elem(neg_idx), beta.elem(neg_idx)) - q.elem(neg_idx);
    neg_idx = arma::find(f_U < 0);
  }

  // Expand the lower bound until f(x_L) <= 0 for all
  arma::uvec pos_idx = arma::find(f_L > 0);
  while (!pos_idx.is_empty()) {
    x_U.elem(pos_idx) = x_L.elem(pos_idx);
    f_U.elem(pos_idx) = f_L.elem(pos_idx); 
    x_L.elem(pos_idx) /= 2;
    f_L.elem(pos_idx) = compound_Poisson_Gamma(x_L.elem(pos_idx), lambda.elem(pos_idx), alpha.elem(pos_idx), beta.elem(pos_idx)) - q.elem(pos_idx);
    pos_idx = arma::find(f_L > 0);
  }

  // Now perform bisection; use an indicator to track active indices
  arma::uvec active = arma::regspace<arma::uvec>(0, n - 1);
  arma::vec x_mid(n, arma::fill::zeros);
  for (int iter = 0; iter < max_iter; ++iter) {

    // Evaluate the central element and check if it is close to the root
    x_mid.elem(active) = 0.5 * (x_L.elem(active) + x_U.elem(active));
    arma::vec f_mid = compound_Poisson_Gamma(x_mid.elem(active), lambda.elem(active), alpha.elem(active), beta.elem(active)) - q.elem(active);
    arma::uvec converged = arma::find(arma::abs(f_mid) < tol);
    if (!converged.is_empty()) {
      active = active.elem(arma::find(arma::abs(f_mid) >= tol));
      if (active.is_empty()) break;
    }

    // Set the middle as lower bound if it is negative
    arma::uvec neg_idx = arma::find(f_mid < 0);
    if (!neg_idx.is_empty()) {
      x_L.elem(active.elem(neg_idx)) = x_mid.elem(neg_idx);
    }

    // Set the middle as upper bound if it is positivd 
    arma::uvec pos_idx = arma::find(f_mid > 0);
    if (!pos_idx.is_empty()) {
      x_U.elem(active.elem(pos_idx)) = x_mid.elem(pos_idx);
    }
  }

  // Mark items where we did not converge as NA and return
  arma::vec result = x_mid;
  if (!active.is_empty()) {
    result.elem(active).fill(R_NaReal);
  }
  return result;
}


// Vectorized Newton-Raphson root-finder for CDF inversion
arma::vec solver_Newton_Raphson(const arma::vec& q, const arma::vec& lambda, const arma::vec& alpha,
                                const arma::vec& beta, const arma::vec& mean, const arma::vec& dispersion,
                                const arma::vec& power) {

  // Set some parameters
  int n = q.n_elem;
  double tol = 1e-8;
  int max_iter = 20;

  // Use the initial point witha woring heurisitic
  arma::vec x0 = mean + 2 * arma::sqrt(dispersion % arma::pow(mean, power)) % (2 * q - 1);
  arma::uvec neg_idx = arma::find(x0 <= 0);
  if (!neg_idx.is_empty()) {
    x0.elem(neg_idx) = mean.elem(neg_idx);
  }
  x0 = arma::clamp(x0, tol, 1 / tol);
  arma::vec x = x0;

  // Keep track of active indices and converged one
  arma::uvec active = arma::regspace<arma::uvec>(0, n - 1);
  arma::vec converged(n, arma::fill::zeros);
  for (int iter = 0; iter < max_iter && !active.is_empty(); ++iter) {
    arma::vec x_active = x.elem(active);

    // For active indices perform a Newton-Raphson iteration (bounding the results)
    arma::vec cdf_x = compound_Poisson_Gamma(x_active, lambda.elem(active), alpha.elem(active), beta.elem(active));
    arma::vec pdf_x = tweedieDensity(x_active, mean.elem(active), dispersion.elem(active), power.elem(active), false);
    arma::vec f = cdf_x - q.elem(active);
    arma::vec f_prime = pdf_x;
    f_prime.elem(arma::find(f_prime == 0)).fill(tol);
    x.elem(active) -= f / f_prime;
    x.elem(active) = arma::clamp(x.elem(active), tol, 1 / std::sqrt(tol));

    // Check convergence for active indices and set converged as non-active
    arma::uvec newly_converged = arma::find(arma::abs(f) < tol);
    if (!newly_converged.is_empty()) {
      converged.elem(active.elem(newly_converged)).ones();
      arma::uvec still_active = arma::find(converged.elem(active) == 0);
      active = active.elem(still_active);
    }
  }

  // After max_iter, set non-converged values using bisection solver and return the root
  if (!active.is_empty()) {
    x.elem(active) = solver_bisection(q.elem(active), lambda.elem(active), alpha.elem(active), beta.elem(active),
                                      x0.elem(active));
  }
  return x;
}


// [[Rcpp::export]]
arma::vec tweedieInvCDF(arma::vec q, arma::vec mean, arma::vec dispersion, arma::vec power) {

  // Make all vectors on the size
  int l = std::max({static_cast<int>(q.n_elem), static_cast<int>(mean.n_elem),
                    static_cast<int>(dispersion.n_elem), static_cast<int>(power.n_elem)});
  q = recycle_to_length(q, l, "q");
  mean = recycle_to_length(mean, l, "mean");
  dispersion = recycle_to_length(dispersion, l, "dispersion");
  power = recycle_to_length(power, l, "power");
  arma::vec invcdf(l, arma::fill::zeros);

  // Get the values of the alternative parametrization
  arma::vec lambda = arma::pow(mean, 2 - power) / (dispersion % (2 - power));
  arma::vec alpha = (2 - power) / (power - 1);
  arma::vec beta = 1 / (dispersion % (power - 1) % arma::pow(mean, power - 1));

  // For stricly positive use teh root-finder
  arma::uvec pos_idx = arma::find(q > arma::exp(-lambda));
  if (!pos_idx.is_empty()) {
    invcdf.elem(pos_idx) = solver_Newton_Raphson(q(pos_idx), lambda(pos_idx), alpha(pos_idx), beta(pos_idx),
                        mean(pos_idx), dispersion(pos_idx), power(pos_idx));

  }
  return invcdf;
}
