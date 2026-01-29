#' Normalize network matrix using double Z-score normalization
#' @param X A numeric matrix (dense or sparse)
#' @return A normalized Matrix object
#' @details Applies row-wise and column-wise Z-score normalization,
#'   combining both with equal weights (1/sqrt(2) each).
#' @noRd
normalizeNetwork <- function(X) {
  nr <- nrow(X)
  nc <- ncol(X)
  
  # Convert to dense matrix for consistent computation
  X_dense <- as.matrix(X)
  
  # Global statistics
  mu0 <- mean(X_dense)
  std0 <- sd(as.vector(X_dense)) * sqrt((nr * nc - 1) / (nr * nc))
  
  # Row-wise Z-scores
  mu1 <- rowMeans(X_dense)
  std1 <- rowSds(X_dense) * sqrt((nc - 1) / nc)
  std1[std1 == 0] <- 1  # Avoid division by zero
  Z1 <- (X_dense - mu1) / std1
  
  # Column-wise Z-scores
  mu2 <- colMeans(X_dense)
  std2 <- colSds(X_dense) * sqrt((nr - 1) / nr)
  std2[std2 == 0] <- 1  # Avoid division by zero
  Z2 <- t((t(X_dense) - mu2) / std2)
  
  # Combined normalization
  normMat <- Z1 / sqrt(2) + Z2 / sqrt(2)
  Z0 <- (X_dense - mu0) / std0
  
  # Handle NA values from zero variance rows/columns
  f1 <- is.na(Z1)
  f2 <- is.na(Z2)
  normMat[f1] <- Z2[f1] / sqrt(2) + Z0[f1] / sqrt(2)
  normMat[f2] <- Z1[f2] / sqrt(2) + Z0[f2] / sqrt(2)
  normMat[f1 & f2] <- 2 * Z0[f1 & f2] / sqrt(2)
  
  # Handle any remaining NA/NaN
  normMat[is.na(normMat)] <- 0
  
  Matrix::Matrix(normMat)
}
