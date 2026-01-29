#' Normalize network matrix using double Z-score normalization
#' @param X A numeric matrix (dense or sparse)
#' @return A normalized Matrix object
#' @details Applies row-wise and column-wise Z-score normalization,
#'   combining both with equal weights (1/sqrt(2) each).
#' @noRd
normalizeNetwork <- function(X) {
  nr <- nrow(X)
  nc <- ncol(X)
  is_sparse <- inherits(X, "sparseMatrix")
  
  # Global statistics
  mu0 <- mean(X)
  std0 <- sd(as.vector(as.matrix(X))) * sqrt((nr * nc - 1) / (nr * nc))
  
  # Row-wise Z-scores
  mu1 <- rowMeans(X)
  std1 <- rowSds(as.matrix(X)) * sqrt((nc - 1) / nc)
  mu1_mat <- matrix(rep(mu1, nc), nrow = nr, ncol = nc)
  std1_mat <- matrix(rep(std1, nc), nrow = nr, ncol = nc)
  Z1 <- (X - mu1_mat) / std1_mat
  
  # Column-wise Z-scores
  mu2 <- colMeans(X)
  std2 <- colSds(as.matrix(X)) * sqrt((nr - 1) / nr)
  mu2_mat <- matrix(rep(mu2, each = nr), nrow = nr, ncol = nc)
  std2_mat <- matrix(rep(std2, each = nr), nrow = nr, ncol = nc)
  Z2 <- (X - mu2_mat) / std2_mat
  
  # Combined normalization
  normMat <- Z1 / sqrt(2) + Z2 / sqrt(2)
  Z0 <- (X - mu0) / std0
  
  # Handle NA values
  if (is_sparse && inherits(normMat, "sparseMatrix")) {
    f1 <- is.na(Z1@x)
    f2 <- is.na(Z2@x)
    Z2_vec <- Z2@x
    Z1_vec <- Z1@x
    Z0_vec <- Z0@x
    normMat@x[f1] <- Z2_vec[f1] / sqrt(2) + Z0_vec[f1] / sqrt(2)
    normMat@x[f2] <- Z1_vec[f2] / sqrt(2) + Z0_vec[f2] / sqrt(2)
    normMat@x[f1 & f2] <- 2 * Z0_vec[f1 & f2] / sqrt(2)
  } else {
    normMat <- as.matrix(normMat)
    Z1 <- as.matrix(Z1)
    Z2 <- as.matrix(Z2)
    Z0 <- as.matrix(Z0)
    f1 <- is.na(Z1)
    f2 <- is.na(Z2)
    normMat[f1] <- Z2[f1] / sqrt(2) + Z0[f1] / sqrt(2)
    normMat[f2] <- Z1[f2] / sqrt(2) + Z0[f2] / sqrt(2)
    normMat[f1 & f2] <- 2 * Z0[f1 & f2] / sqrt(2)
    normMat <- Matrix::Matrix(normMat)
  }
  
  return(normMat)
}
