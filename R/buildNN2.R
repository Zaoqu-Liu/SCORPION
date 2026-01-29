#' Build k-nearest neighbor graph using nn2
#' @param X A coordinate matrix (samples x features)
#' @param k Number of nearest neighbors
#' @param mode Graph mode for igraph
#' @return A list containing the kNN graph
#' @importFrom RANN nn2
#' @importFrom igraph graph_from_adj_list simplify
#' @noRd
buildNN2 <- function(X, k = min(5, ncol(X)), mode = "all") {
  nn2.res <- RANN::nn2(data = X, k = k)$nn.idx
  adj.knn <- split(nn2.res, rep(1:nrow(nn2.res), times = ncol(nn2.res)))
  
  graph.knn <- igraph::graph_from_adj_list(adj.knn, duplicate = FALSE, mode = mode)
  graph.knn <- igraph::simplify(graph.knn, remove.multiple = TRUE)
  igraph::E(graph.knn)$weight <- 1
  
  list(graph.knn = graph.knn)
}
