#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List eureka(int order_max, vec r, vec g, mat coefs, vec var, vec a, double threshold)  {
  
  int l, l1, l2, i, j, k;
  double v, d, q, hold;
  List ans = List::create(_["vars"] = var, _["coefs"] = coefs, _["order"] = order_max);
  
  v = r(0);
  d = r(1);
  a(0) = 1.0;
  coefs(0,0) = g(1)/v;
  
  if(abs(coefs(0,0)) <= threshold) {
    ans["vars"] = var;
    ans["coefs"] = coefs;
    ans["order"] = 0;
    return ans;
  }
  
  q = coefs(0,0) * r(1);
  var(0) = (1 - coefs(0,0)*coefs(0,0)) * r(0);
  
  if(order_max == 1) {
    ans["vars"] = var;
    ans["coefs"] = coefs;
    ans["order"] = 1;
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
    
    if(abs(coefs(l-1, l-1)) <= threshold) {
      ans["vars"] = var;
      ans["coefs"] = coefs;
      ans["order"] = l-1;
      return ans;
    }
    
    // since coefficienets will most likely be in decreasing order, coefs[l-1,l-1] will be the smallest.
    // Therefore checking the first occurence of where this is smaller than the threshold would be enough.
    
    for(j = 1; j <= l-1; j++) {
      coefs(l-1, j-1) = coefs(l-2, j-1) + coefs(l-1, l-1) * a(l-j);
    }
    
    var(l-1) = var(l-2) * (1 - coefs(l-1,l-1)*coefs(l-1,l-1));
    
    if(l == order_max) {
      ans["vars"] = var;
      ans["coefs"] = coefs;
      ans["order"] = order_max;
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
  ans["order"] = order_max;
  return ans;
}


List ar_yw(int order_max, vec r, vec g, mat coefs, vec var, vec a, double threshold, uword n)  {
  List ans = eureka(order_max, r, g, coefs, var, a, threshold);
  
  vec vars = ans["vars"];
  mat coef_mat = ans["coefs"];
  int order = ans["order"];

  vec coef_vec(order, fill::zeros);
  double var_pred = 0;
  
  if(order> 0)  {
    coef_vec = coef_mat(order-1, span(0, order-1)).t();
    var_pred = vars(order-1);
  }
  else{
    var_pred = r(0);
  }
  
  var_pred = var_pred*n/(n-(order+1));
  
  List ret = List::create(_["coefs"] = coef_vec, _["vars"] = var_pred, _["order"] = order);
  return(ret);
  
}

List arp_approx(vec xacf, int max_order, uword n, double threshold) {
  List ar_fit = ar_yw(max_order, xacf, xacf, zeros<mat>(max_order, max_order),
                      zeros<vec>(max_order), zeros<vec>(max_order), threshold, n);
  
  vec coefs = ar_fit["coefs"];
  double vars = ar_fit["vars"];
  int order = ar_fit["order"];
  
  double spec = vars/pow(1 - sum(coefs),2.0);
  double foo = 0, Gamma = 0;
  
  if(order != 0)  {
    foo = 0;
    
    for(int i=1; i<= order; i++) {
      for(int k=1; k<=i; k++) {
        foo = foo + coefs(i-1) * k * xacf(abs(k-i));
      }
    }
    
    Gamma = 2*(foo + (spec-xacf(0))/2 * dot(regspace(1, order), coefs)) / (1 - sum(coefs));
  }
  else
    Gamma = 0;
  
  List ans = List::create(_["Gamma"] = Gamma, _["Sigma"] = spec);
  return(ans);
}


// [[Rcpp::export]]
double batchsize_cpp(uword n, int p, mat xacf_mat, int max_order, String method, double threshold = 0.01) {

  mat ar_fit(2,p, fill::none);
  double num_sum = 0, denom_sum = 0, Gamma = 0, Sigma = 0;
  List arp_output = List::create(_["Gamma"] = 0, _["Sigma"] = 0);
  
  
  for(int i = 0; i<p; i++)  {

    arp_output = arp_approx(xacf_mat.col(i), max_order, n, threshold);
    Gamma = arp_output["Gamma"];
    Sigma = arp_output["Sigma"];
    ar_fit(0,i) = pow(Gamma,2.0);
    ar_fit(1,i) = pow(Sigma,2.0);
    num_sum += ar_fit(0,i);
    denom_sum += ar_fit(1,i);
  }
  
  double coeff = pow(num_sum/denom_sum,1.0/3);
  double b_const = (3.0/2*n)*(method == "obm" || method == "bartlett" || method == "tukey")
    + (n)*(method == "bm");
  double b = pow(b_const,1.0/3) * coeff;
  if(b < 1)
    b = 1;
  b = floor(b);
  return(coeff);
}

