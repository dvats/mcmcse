#include <RcppArmadillo.h>
#include <math.h> 

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat mobmC(const arma::mat& chain, double b)
{
  int n = chain.n_rows;
  int p = chain.n_cols;
  int a = n-b+1;

  // y_mean will store the overall mean vector
  vec y_mean(p);


  // out is the output matrix that returns the Sigma_hat
  // block_means stores the mean in a block of size b, overall a blocks
  // mean_mat is the matrix will each column being y_mean. This matrix is made for easy calculation
  mat out(p,p);
  mat block_means(a,p);
  mat mean_mat(a,p);

  block_means.zeros();
  out.zeros();
  y_mean.zeros();
  mean_mat.zeros();

  // idx will be used to find block_means. Helps bypass a double loop
  IntegerVector foo = seq_len(a);
  uvec idx = as<uvec>(foo) - 1;
  // idx = idx*b;

 // Find Block means
  for(int i = 0; i < b; ++i)
  {
    block_means += chain.rows(idx);
    idx += 1;
  }
  block_means = block_means/b;
   

 // Find the overall mean vector
  for(int i = 0; i < n; ++i)
  {
    y_mean += trans(chain.rows(i,i));
  }
  y_mean = y_mean/n;

// make the matrix of mean vector
  for(int i = 0; i < a; i++)
  {
    mean_mat.row(i) = trans(y_mean);
  }

  out += trans(block_means - mean_mat)*(block_means - mean_mat);
  return(out*b/n);
}