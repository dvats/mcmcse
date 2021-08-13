#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


uword counting_obm(vec var_vector, double var_number) {
  uword count = 0;
  for(uword i = 0; i<var_vector.size(); i++) {
    if(var_vector(i) <=  var_number)
      count++;
  }
  return(count);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mcseqobm(const arma::vec& x, double b, double xi_hat) {
  uword n = x.n_elem;
  double a = n-b+1;
  vec y(a);
  for(uword k=1; k<=a; k++)  {
    y(k-1) = counting_obm( x(span((k-1),(k + b - 2))), xi_hat);
  }
  y = y/b;
  double mu_hat = mean(y);
  double var_hat = n * b * sum( pow(y - mu_hat, 2.0)) / (a - 1.0) / a;
  return(var_hat);
}

