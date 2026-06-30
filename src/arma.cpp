#include <Rcpp.h>
using namespace Rcpp;

// The dynamic estimation of the parameters of a Gamma distribution
// [[Rcpp::export]]
List armaDynamic(NumericVector y, double phi, double theta, double co) {
  int n = y.size();
  NumericVector m(n);
  NumericVector v(n);
  
  double f = var(y);
  double K = phi+theta;
  
  m[0] = 0;
  v[0] = y[1];
  if (n >= 2) {
    for (int t = 1; t < n; t++) {
      m[t] = phi * m[t-1] + K * v[t-1];
      v[t] = y[t] - co - m[t];
    }
  }
  
  return List::create(_["m"] = m, _["v"] = v);
}