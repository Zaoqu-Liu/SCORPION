#' Convert edge list to adjacency matrix
#' @param df A data.frame with 3 columns: source, target, weight
#' @return A data.frame representing the adjacency matrix
#' @noRd
spreadNet <- function(df) {
  df[, 3] <- as.numeric(df[, 3])
  row_names <- unique(df[, 1])
  col_names <- unique(df[, 2])
  
  # Use matrix indexing instead of loop
  spread.df <- matrix(0, nrow = length(row_names), ncol = length(col_names),
                      dimnames = list(row_names, col_names))
  idx <- cbind(match(df[, 1], row_names), match(df[, 2], col_names))
  spread.df[idx] <- df[, 3]
  
  as.data.frame(spread.df)
}
