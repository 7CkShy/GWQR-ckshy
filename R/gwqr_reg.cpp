// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>


using namespace Rcpp;
using namespace arma;
using namespace stats;


// [[Rcpp::export]]
mat gwqr_uv(mat dp, mat X){
  int n = dp.n_rows;
  int p = X.n_cols;
  arma::mat U,V;
  U.zeros(n,n);
  //arma::mat V;
  V.zeros(n,n);
  for (int j=0; j < n; j++){
    for (int i=0; i < n; i++){
      U(i,i) = abs(dp(i, 0) - dp(j, 0));
      V(i,i) = abs(dp(i, 1) - dp(j, 1));
    }
  }
  mat XX = X;
  X.insert_cols(p, U*XX);
  X.insert_cols(p+p, V*XX);
  return X;
}



// [[Rcpp::export]]
mat gwqr_reg(mat dp, mat X, vec Y){
  int n = X.n_rows, p = X.n_cols;
  arma::mat A;
  for (int i = 0; i < n; i++){
    A(i,i) = 1;
  }
  return A;
}