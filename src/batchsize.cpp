#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
// using namespace arma;

List eureka(int order_max, const arma::vec& r, const arma::vec& g, arma::mat coefs, arma::vec var, arma::vec a, double threshold)  {
// Using the levinson algorithm to calculate AR coefficients and variance.
// We have added a check which compares the coefficient to a threshold. If at any time a coefficient falls 
// below the threshold, we don't compute coefficients beyond that order (i.e. we don't always complete the 
// loop from 1 to order.max) and return. 

  int l, i, j, k;
  double v, d, q, hold, l1, l2;
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
  var(0) = (1.0 - coefs(0,0)*coefs(0,0)) * r(0);

  if(order_max == 1) {
    ans["vars"] = var;
    ans["coefs"] = coefs;
    ans["order"] = 1;
    return ans;
  }

  for(l=2; l<= order_max; l++) {
    a(l-1) = -d/v;

    if(l > 2) {
      l1 = (l-2.0)/2.0;
      l2 = l1+1.0;

      for(j = 2; j<= l2; j++) {
        hold = a(j-1);
        k = l-j+1;
        a(j-1) = a(j-1) + a(l-1) * a(k-1);
        a(k-1) = a(k-1) + a(l-1) * hold;
      }

      if((2.0*l1) != (l-2))
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

    var(l-1) = var(l-2) * (1.0 - coefs(l-1,l-1)*coefs(l-1,l-1));

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


List ar_yw(int order_max, const arma::vec& r, const arma::vec& g, arma::mat coefs, arma::vec var, arma::vec a, double threshold, arma::uword n)  {
// Wrapper function to calculate AR coefficients and variance from the output of Levinson algorithm (eureka).
  List ans = eureka(order_max, r, g, coefs, var, a, threshold);

  arma::vec vars = ans["vars"];
  arma::mat coef_mat = ans["coefs"];
  int order = ans["order"];
  arma::vec coef_vec(order, arma::fill::zeros);
  double var_pred = 0.0;

  if(order> 0)  {
    coef_vec = coef_mat(order-1, arma::span(0, order-1)).t();
    var_pred = vars(order-1);
  }
  else{
    var_pred = r(0);
  }

  var_pred = var_pred*n/(n-(order+1.0));

  List ret = List::create(_["coefs"] = coef_vec, _["vars"] = var_pred, _["order"] = order);
  return(ret);

}

List arp_approx(const arma::vec& xacf, int max_order, arma::uword n, double threshold) {
// Calculate Gamma and Sigma from AR approximation
  List ar_fit = ar_yw(max_order, xacf, xacf, arma::zeros<arma::mat>(max_order, max_order),
                      arma::zeros<arma::vec>(max_order), arma::zeros<arma::vec>(max_order), threshold, n);

  arma::vec coefs = ar_fit["coefs"];
  double vars = ar_fit["vars"];
  int order = ar_fit["order"];

  double spec = vars/pow(1.0 - sum(coefs),2.0);
  double foo = 0.0, Gamma = 0.0;

  if(order != 0)  {
    foo = 0.0;

    for(int i=1; i<= order; i++) {
      for(int k=1; k<=i; k++) {
        foo = foo + coefs(i-1) * k * xacf(abs(k-i));   // convert k to floor?
      }
    }

    Gamma = 2.0*(foo + (spec-xacf(0))/2.0 * dot(arma::regspace(1, order), coefs)) / (1.0 - sum(coefs));
  }
  else
    Gamma = 0.0;

  List ans = List::create(_["Gamma"] = Gamma, _["Sigma"] = spec);
  return(ans);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double batchsize_cpp(arma::uword n, int p, const arma::mat& xacf_mat, int max_order, String method, double threshold = 0.01) {
// Calculates batchsize from Gamma and Sigma from AR approximation. Uses the "optimal" batchsize parametric 
// method  from of Liu et.al.
  arma::mat ar_fit(2,p, arma::fill::none);
  double num_sum = 0.0, denom_sum = 0.0, Gamma = 0.0, Sigma = 0.0;
  List arp_output = List::create(_["Gamma"] = Gamma, _["Sigma"] = Sigma);

  for(int i = 0; i<p; i++)  
  {
    arp_output = arp_approx(xacf_mat.col(i), max_order, n, threshold);
    Gamma = arp_output["Gamma"];
    Sigma = arp_output["Sigma"];
    ar_fit(0,i) = pow(Gamma,2.0);
    ar_fit(1,i) = pow(Sigma,2.0);
    num_sum += ar_fit(0,i);
    denom_sum += ar_fit(1,i);
  }

  double coeff = pow(num_sum/denom_sum, (1.0/3.0) );
  double b_const = ( (3.0/2.0)*n)*(method == "obm" || method == "bartlett" || method == "tukey" || method == "sub")
                  + (n)*(method == "bm");
  double b = pow(b_const,1.0/3.0) * coeff;
  if(b < 1.0)
    b = 1.0;
  // b = floor(b);
  return(b);
}

