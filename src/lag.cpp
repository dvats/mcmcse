#include <RcppArmadillo.h>
#include <math.h> 

// using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec lag(arma::vec s, int n, double b, String method)
{
  arma::vec w(n);
  w.zeros();
  if(method == "bartlett")
  {
    w.head(b) = 1 - s/b;
  }
  else if(method == "tukey")
  {
    w.head(b) = (1 + cos(M_PI * s/b))/2 ;
  }
  else
  {
    stop("Invalid method. Only bartlett and tukey allowed");
  }
  return(w);
}
