#' Column-wise standard deviation
#' @param X A numeric matrix
#' @param na.rm Logical, whether to remove NA values
#' @return A numeric vector of standard deviations
#' @noRd
colSds <- function(X, na.rm = FALSE) {
  if (na.rm) {
    n <- colSums(!is.na(X))
    means <- colSums(X, na.rm = TRUE) / n
    sq_diff <- sweep(X, 2, means, "-")^2
    sqrt(colSums(sq_diff, na.rm = TRUE) / (n - 1))
  } else {
    n <- nrow(X)
    means <- colMeans(X)
    sq_diff <- sweep(X, 2, means, "-")^2
    sqrt(colSums(sq_diff) / (n - 1))
  }
}
