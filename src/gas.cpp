#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

template <typename ScoreFn, typename InvLinkFn>
List gasFilterCore(NumericVector y, double psi0, double phi, double rho,
                   double xi0, double k, int period, ScoreFn score_fn, InvLinkFn inv_link_fn) {
  int n = y.size();
  if (n == 0) stop("y is empty");
  if (period < 1) stop("period must be >= 1");

  NumericVector f(n);
  NumericVector xi(period);
  for (int j = 0; j < period; ++j) xi[j] = xi0;

  double psi = psi0;

  for (int t = 0; t < n; ++t) {
    double score = 0.0;
    if (t > 0) score = score_fn(y[t - 1], f[t - 1]);

    psi = phi * psi + rho * score;

    if (period > 1) {
      int past_season = (t - 1) % period;
      xi[past_season] += k * score;

      double adj = k * score / (period - 1.0);
      for (int j = 0; j < period; ++j) {
        if (j != past_season) xi[j] -= adj;
      }

      int current_season = t % period;
      f[t] = inv_link_fn(psi + xi[current_season]);
    } else {
      f[t] = inv_link_fn(psi);
    }
  }

  return List::create(
    _["f"] = f,
    _["last_psi"] = psi,
    _["last_xi"] = xi
  );
}


// [[Rcpp::export]]
List gasFilterPois(NumericVector y, double psi0, double phi, double rho,
                   double xi0, double k, int period) {
  auto score_pois = [](double y, double f) { return y - f; };
  auto inv_link_log = [](double x) { return std::exp(x); };
  return gasFilterCore(y, psi0, phi, rho, xi0, k, period, score_pois, inv_link_log);
}

// [[Rcpp::export]]
List gasFilterNbinom(NumericVector y, double psi0, double phi, double rho,
                     double xi0, double k, int period, double alpha) {
  if (!std::isfinite(alpha) || alpha <= 0.0) stop("alpha must be > 0");
  auto score_nbinom = [alpha](double y, double f) { return (y - f) / (1.0 + f / alpha);};
  auto inv_link_log = [](double x) { return std::exp(x); };
  return gasFilterCore(y, psi0, phi, rho, xi0, k, period, score_nbinom, inv_link_log);
}

// [[Rcpp::export]]
List gasFilterBern(NumericVector y, double psi0, double phi, double rho,
                   double xi0, double k, int period) {
    auto score_bern = [](double y, double f) { return (1.0 - y) - f; };
    auto inv_link_logit = [](double x) { return 1.0 / (1.0 + std::exp(-x)); };
    return gasFilterCore(y, psi0, phi, rho, xi0, k, period, score_bern, inv_link_logit);
}






