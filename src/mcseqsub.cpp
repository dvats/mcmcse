#include <iostream>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// uword counting(vec var_vector, double var_number) {
//   uword count = 0;
//   for(uword i = 0; i<var_vector.size(); i++) {
//     if(var_vector(i) <=  var_number)
//       count++;
//   }
//   return(count);
// }

// [[Rcpp::export]]
double mcseqsub(const arma::vec& x, double b, double q, Function f) {
  uword n = x.n_elem;
  double a = n-b+1;
  vec P = {q};
  vec y(a);
  for(uword k=1; k<=a; k++)  {
    y(k-1) = *(REAL(f(x(span((k-1),(k + b - 2))), q)));
  }
  double mu_hat = mean(y);
  double var_hat = n * b * sum( pow(y - mu_hat, 2.0)) / (a - 1.0) / a;
  return(var_hat);
}

