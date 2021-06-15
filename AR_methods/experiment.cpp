#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec fun(mat chain) {
  Environment stats("package:stats");
  Function acf = stats["acf"];
  Function na_pass = stats["na.pass"];
  List result =  acf(chain.col(1), 5, "covariance", false, na_pass);
  vec c = result["acf"];
  //Rcout << c;
  // Rcpp::Environment base("package:stats"); 
  // Function acf = base["acf"];
  // RObject v = acf(_["x"] = x, _["type"] = "covariance", _["lag.max"] = 5, _["plot"] = FALSE, _["demean"]=TRUE);
  return(result["acf"]);
}