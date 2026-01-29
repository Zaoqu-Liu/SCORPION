#' Row-wise standard deviation
#' @param X A numeric matrix
#' @param na.rm Logical, whether to remove NA values
#' @return A numeric vector of standard deviations
#' @noRd
rowSds <- function(X, na.rm = FALSE) {
  if (na.rm) {
    n <- rowSums(!is.na(X))
    means <- rowSums(X, na.rm = TRUE) / n
    sq_diff <- sweep(X, 1, means, "-")^2
    sqrt(rowSums(sq_diff, na.rm = TRUE) / (n - 1))
  } else {
    n <- ncol(X)
    means <- rowMeans(X)
    sq_diff <- sweep(X, 1, means, "-")^2
    sqrt(rowSums(sq_diff) / (n - 1))
  }
}
