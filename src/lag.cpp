#include <RcppArmadillo.h>
#include <math.h> 

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec lag(vec s, int n, double b, String method)
{
  vec w(n);
  w.zeros();
  if(method == "bartlett")
  {
    w.head(b) = 1 - s/b;
  }
  else if(method == "tukey")
  {
    w.head(b) = (1 + cos(PI * s/b))/2 ;
  }
  else
  {
    stop("Invalid method. Only bartlett and tukey allowed");
  }
  return(w);
}
