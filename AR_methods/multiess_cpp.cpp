#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double multiESS_cpp(mat chain, mat covmat)  {
  uword n = size(chain)(0);
  int p = size(chain)(1);
  double det_var=0, det_cov=0;
  mat var_mat = cov(chain);
  vec eigval_var = pow(eig_sym(var_mat), 1.0/p);
  vec eigval_cov = pow(eig_sym(covmat), 1.0/p);
  det_var = prod(eigval_var);
  det_cov = prod(eigval_cov);
  double ess = n * (det_var/det_cov);
  return(ess);
}