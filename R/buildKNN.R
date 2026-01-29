#' Build k-nearest neighbor graph
#' @param X Input data (coordinates or distance matrix)
#' @param k Number of nearest neighbors
#' @param from Input type: "dist" or "coordinates"
#' @param use.nn2 Whether to use RANN::nn2 for speed
#' @param return_neighbors_order Whether to return neighbor indices
#' @param dist_method Distance metric
#' @param cor_method Correlation method if dist_method is "cor"
#' @param p Minkowski p parameter
#' @param directed Whether to build directed graph
#' @return A list containing the kNN graph
#' @noRd
buildKNN <- function(X,
                     k = 5,
                     from = c("dist", "coordinates"),
                     use.nn2 = TRUE,
                     return_neighbors_order = FALSE,
                     dist_method = "euclidean",
                     cor_method = "pearson",
                     p = 2,
                     directed = FALSE) {
  av.methods <- c("dist", "coordinates")
  method <- pmatch(from[1], av.methods)
  if (is.na(method)) {
    stop(paste("Unknown method:", from, 
               "\nAvailable methods:", paste(av.methods, collapse = ", ")))
  }
  
  if (method == 2) {
    # From coordinates
    if (use.nn2) {
      if (dist_method != "euclidean") {
        stop("nn2 only supports euclidean distance. Set use.nn2 = FALSE for other distances.")
      }
      mode <- ifelse(directed, "out", "all")
      return(buildNN2(X = X, k = k, mode = mode))
    } else {
      av.dist <- c("cor", "euclidean", "maximum", "manhattan", 
                   "canberra", "binary", "minkowski")
      dist_method_ <- pmatch(dist_method, av.dist)
      if (is.na(dist_method_)) {
        stop(paste("Unknown distance method:", dist_method,
                   "\nAvailable:", paste(av.dist, collapse = ", ")))
      }
      if (dist_method_ == 1) {
        av.cor_methods <- c("pearson", "kendall", "spearman")
        cor_method_ <- pmatch(cor_method, av.cor_methods)
        if (is.na(cor_method_)) {
          stop(paste("Unknown correlation method:", cor_method,
                     "\nAvailable:", paste(av.cor_methods, collapse = ", ")))
        }
        X <- stats::as.dist(as.matrix(1 - stats::cor(t(X), method = cor_method)))
      } else {
        X <- stats::dist(X, method = dist_method)
      }
    }
  } else {
    if (use.nn2) {
      stop("nn2 requires coordinates, not distance. Set from = 'coordinates' or use.nn2 = FALSE")
    }
    return(buildKNND(D = X, k = k, return_neighbors_order = return_neighbors_order))
  }
  
  buildKNND(D = X, k = k, return_neighbors_order = return_neighbors_order)
}
