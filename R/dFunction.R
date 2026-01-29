#' Compute non-negative correlation matrix
#' @param X First matrix
#' @param Y Second matrix
#' @return Correlation matrix with negative values set to zero
#' @noRd
dFunction <- function(X, Y) {
  A <- cor(X, Y)
  A[A < 0] <- 0
  A
}
