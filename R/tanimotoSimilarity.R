#' Compute Tanimoto similarity between two matrices
#' @description Computes a modified Tanimoto coefficient between row vectors of X
#'   and column vectors of Y. Used in PANDA message passing.
#'   Uses RcppArmadillo for optimized BLAS-accelerated computation.
#' @param X A numeric matrix (rows: features, cols: samples)
#' @param Y A numeric matrix (rows: samples, cols: features)
#' @return A similarity matrix of dimension nrow(X) x ncol(Y)
#' @useDynLib SCORPION, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @noRd
tanimoto <- function(X, Y) {
  # Use C++ implementation for speed
  tanimoto_cpp(as.matrix(X), as.matrix(Y))
}
