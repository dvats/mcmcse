#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
uword counting(vec var_vector, double var_number) {
  uword count = 0;
  for(uword i = 0; i<var_vector.size(); i++) {
    if(var_vector(i) <=  var_number)
      count++;
  }
  return(count);
}

// [[Rcpp::export]]
double mcseqbm(const arma::vec& x, double b, double xi_hat) {
  uword n = x.n_elem;
  double a = floor(n/b);
  vec y(a);
  for(uword k=1; k<=a; k++)  {
    y(k-1) = counting( x(span(((k - 1) * b ),(k * b - 1))), xi_hat);
  }
  y = y/b;
  double mu_hat = mean(y);
  double var_hat = b * sum( pow(y - mu_hat, 2.0)) / (a - 1.0);
  return(var_hat);
}

