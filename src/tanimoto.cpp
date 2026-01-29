// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Fast Tanimoto similarity computation using RcppArmadillo
//' @param X First matrix (rows: features, cols: samples)
//' @param Y Second matrix (rows: samples, cols: features)
//' @return Tanimoto similarity matrix
//' @noRd
// [[Rcpp::export]]
arma::mat tanimoto_cpp(const arma::mat& X, const arma::mat& Y) {
  int nr = X.n_rows;
  int nc = Y.n_cols;
  
  // Compute dot product using optimized BLAS
  arma::mat Amat = X * Y;
  
  // Column-wise sum of squares for Y
  arma::rowvec Bvec = sum(Y % Y, 0);  // 1 x nc
  
  // Row-wise sum of squares for X
  arma::colvec Cvec = sum(X % X, 1);  // nr x 1
  
  // Expand to matrices
  arma::mat Bmat = repmat(Bvec, nr, 1);  // nr x nc
  arma::mat Cmat = repmat(Cvec, 1, nc);  // nr x nc
  
  // Tanimoto coefficient: A / sqrt(B + C - |A|)
  arma::mat den = Bmat + Cmat - abs(Amat);
  
  // Avoid division by zero
  den.elem(find(den <= 0)).fill(datum::eps);
  
  Amat = Amat / sqrt(den);
  
  return Amat;
}

//' Fast matrix normalization
//' @param X Input matrix
//' @return Normalized matrix with mean 0 and unit L2 norm per row
//' @noRd
// [[Rcpp::export]]
arma::mat normalize_rows_cpp(const arma::mat& X) {
  arma::mat result = X;
  
  // Subtract row means
  arma::colvec row_means = mean(X, 1);
  result.each_col() -= row_means;
  
  // Normalize by L2 norm
  arma::colvec row_norms = sqrt(sum(result % result, 1));
  
  // Avoid division by zero
  row_norms.elem(find(row_norms == 0)).fill(1.0);
  
  result.each_col() /= row_norms;
  
  return result;
}

//' Update diagonal elements during PANDA iteration
//' @param X Input matrix
//' @param n Dimension
//' @param alpha Learning rate
//' @param step Current iteration step
//' @return Matrix with updated diagonal
//' @noRd
// [[Rcpp::export]]
arma::mat update_diagonal_cpp(arma::mat X, int n, double alpha, int step) {
  double val = -static_cast<double>(step * 2 + n) / (step * 2 + 2);
  X.diag().fill(val);
  return X;
}
