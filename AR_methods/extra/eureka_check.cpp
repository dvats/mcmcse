#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Rcpp copy of eureka.f to check that it's answer matches with ar.yw()

// [[Rcpp::export]]
List eureka_check(int order_max, vec r, vec g, mat coefs, vec var, vec a)  {
  
  int l, l1, l2, i, j, k;
  double v, d, q, hold;
  List ans = List::create(_["vars"] = var, _["coefs"] = coefs);
  
  v = r(0);
  d = r(1);
  a(0) = 1.0;
  coefs(0,0) = g(1)/v;
  q = coefs(0,0) * r(1);
  var(0) = (1 - coefs(0,0)*coefs(0,0)) * r(0);
  
  if(order_max == 1) {
    ans["vars"] = var;
    ans["coefs"] = coefs;
    return ans;
  }
  
  for(l=2; l<= order_max; l++) {
    a(l-1) = -d/v;
    
    if(l > 2) {
      l1 = (l-2)/2;
      l2 = l1+1;
      
      for(j = 2; j<= l2; j++) {
        hold = a(j-1);
        k = l-j+1;
        a(j-1) = a(j-1) + a(l-1) * a(k-1);
        a(k-1) = a(k-1) + a(l-1) * hold;
      }
      
      if((2*l1) != (l-2))
        a(l2) = a(l2) * (1.0 + a(l-1));
    }
    
    v = v + a(l-1) * d;
    coefs(l-1, l-1) = (g(l) - q)/v;
    
    for(j = 1; j <= l-1; j++) {
      coefs(l-1, j-1) = coefs(l-2, j-1) + coefs(l-1, l-1) * a(l-j);
    }
    
    var(l-1) = var(l-2) * (1 - coefs(l-1,l-1)*coefs(l-1,l-1));
    
    if(l == order_max) {
      ans["vars"] = var;
      ans["coefs"] = coefs;
      return ans;
    }
    
    d = 0.0;
    q = 0.0;
    
    for(i=1; i<=l; i++) {
      k = l-i+2;
      d = d + a(i-1) * r(k-1);
      q = q + coefs(l-1, i-1) * r(k-1);
    }
    
  }
  
  ans["vars"] = var;
  ans["coefs"] = coefs;
  return ans;
}

