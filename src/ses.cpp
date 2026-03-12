#include <Rcpp.h>
using namespace Rcpp;

// A damped simple exponential smoothing with learnable parmaeters.
// [[Rcpp::export]]
NumericVector dampedSES(NumericVector y, double mu0, double alpha, double phi) {
  int n = y.size();
  if (n <= 1) {
    Rcpp::stop("y is too short: at least 2 observations required");
  }
  if (alpha < 0 || phi < 0) {
    Rcpp::stop("Negative smoothing parameters");
  }

  while (alpha + phi > 1) {
    double param_sum = alpha + phi;
    alpha = alpha / param_sum;
    phi = phi / param_sum;
  }
  double complement = 1 - alpha - phi;

  NumericVector mu(n);
  double ymean = mean(y);

  mu[0] = mu0;
  for (int t = 1; t < n; t++) {
    mu[t] = alpha * y[t - 1] + phi * ymean + complement * mu[t - 1];
  }

  return mu;
}
