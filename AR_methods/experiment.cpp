#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double multiESS_cpp(arma::mat chain, arma::mat covmat)  {
  uword n = size(chain)(0);
  int p = size(chain)(1);
  double log_det_var=0, log_det_cov=0;
  mat var_mat = cov(chain);
  vec log_eigval_var = log(eig_sym(var_mat));
  vec log_eigval_cov = log(eig_sym(covmat));
  log_det_var = sum(log_eigval_var);
  log_det_cov = sum(log_eigval_cov);
  double ess = n * exp((log_det_var - log_det_cov) * (1.0/p));
  return(ess);
}