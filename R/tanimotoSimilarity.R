#' Compute Tanimoto similarity between two matrices
#' @param X A numeric matrix (rows: features, cols: samples)
#' @param Y A numeric matrix (rows: samples, cols: features)
#' @return A similarity matrix
#' @noRd
tanimoto <- function(X, Y) {
  nc <- ncol(Y)
  nr <- nrow(X)
  dm <- c(nr, nc)
  
  # Compute dot product

  Amat <- X %*% Y
  
  # Column-wise sum of squares
  Bmat <- colSums(Y * Y)
  Bmat <- rep(Bmat, each = nr)
  dim(Bmat) <- dm
  
  # Row-wise sum of squares
  Cmat <- rowSums(X * X)
  Cmat <- rep(Cmat, nc)
  dim(Cmat) <- dm
  
  # Tanimoto coefficient

  den <- Bmat + Cmat - abs(Amat)
  Amat <- Amat / sqrt(den)
  
  return(Amat)
}
