#' Compute Tanimoto similarity between two matrices
#' @description Computes a modified Tanimoto coefficient between row vectors of X
#'   and column vectors of Y. Used in PANDA message passing.
#' @param X A numeric matrix (rows: features, cols: samples)
#' @param Y A numeric matrix (rows: samples, cols: features)
#' @return A similarity matrix of dimension nrow(X) x ncol(Y)
#' @noRd
tanimoto <- function(X, Y) {
  nc <- ncol(Y)
  nr <- nrow(X)
  dm <- c(nr, nc)
  
  # Compute dot product
  Amat <- X %*% Y
  
  # Column-wise sum of squares for Y
  Bmat <- colSums(Y * Y)
  Bmat <- rep(Bmat, each = nr)
  dim(Bmat) <- dm
  
  # Row-wise sum of squares for X
  Cmat <- rowSums(X * X)
  Cmat <- rep(Cmat, nc)
  dim(Cmat) <- dm
  
  # Tanimoto coefficient: A / sqrt(B + C - |A|)
  den <- Bmat + Cmat - abs(Amat)
  # Avoid division by zero or negative values
  den[den <= 0] <- .Machine$double.eps
  Amat <- Amat / sqrt(den)
  
  return(Amat)
}
