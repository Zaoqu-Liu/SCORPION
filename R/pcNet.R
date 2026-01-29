#' Principal Component Regression Network
#' @description Constructs a gene co-expression network using principal component regression
#' @param X A gene expression matrix (genes x cells)
#' @param nComp Number of principal components to use
#' @param scaleScores Whether to scale output to [-1, 1]
#' @param symmetric Whether to make the output symmetric
#' @param q Quantile threshold for filtering
#' @param verbose Whether to show progress
#' @param nCores Number of cores for parallel computation
#' @return A sparse Matrix of regression coefficients
#' @importFrom pbapply pbsapply
#' @noRd
pcNet <- function(X,
                  nComp = 3,
                  scaleScores = TRUE,
                  symmetric = FALSE,
                  q = 0,
                  verbose = FALSE,
                  nCores = 1) {
  # Input validation
  if (!all(Matrix::rowSums(X) > 0)) {
    stop("Quality control has not been applied over the matrix.")
  }
  xClass <- class(X)[[1]]
  if (!xClass %in% c("matrix", "dgCMatrix")) {
    stop("Input should be a matrix with cells as columns and genes as rows")
  }
  if (nComp < 2 || nComp >= nrow(X)) {
    stop("nComp should be >= 2 and < total number of genes")
  }
  
  gNames <- rownames(X)
  X <- scale(Matrix::t(X))
  n <- ncol(X)
  
  # Principal component regression for each gene
  pcCoefficients <- function(K) {
    y <- X[, K]
    Xi <- X[, -K]
    coeff <- irlba::irlba(Xi, nComp)$v
    score <- Xi %*% coeff
    score_norms <- sqrt(colSums(score^2))
    score <- sweep(score, 2, score_norms^2, "/")
    Beta <- colSums(y * score)
    coeff %*% Beta
  }
  
  # Initialize output matrix
  A <- 1 - diag(n)
  
  # Compute coefficients with optional parallelization
  if (nCores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(nCores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("X", "nComp"), envir = environment())
      parallel::clusterEvalQ(cl, library(irlba))
      B <- if (verbose) {
        pbapply::pbsapply(seq_len(n), pcCoefficients, cl = cl)
      } else {
        parallel::parSapply(cl, seq_len(n), pcCoefficients)
      }
    } else {
      B <- if (verbose) {
        pbapply::pbsapply(seq_len(n), pcCoefficients, cl = nCores)
      } else {
        simplify2array(parallel::mclapply(seq_len(n), pcCoefficients, mc.cores = nCores))
      }
    }
  } else {
    B <- if (verbose) {
      pbapply::pbsapply(seq_len(n), pcCoefficients)
    } else {
      sapply(seq_len(n), pcCoefficients)
    }
  }
  
  # Fill output matrix
  B <- t(B)
  diag_mask <- diag(n) == 0
  A[diag_mask] <- as.vector(B)
  
  # Post-processing
  if (isTRUE(symmetric)) {
    A <- (A + t(A)) / 2
  }
  
  absA <- abs(A)
  if (isTRUE(scaleScores)) {
    A <- A / max(absA)
  }
  
  A[absA < stats::quantile(absA, q)] <- 0
  diag(A) <- 0
  colnames(A) <- rownames(A) <- gNames
  
  Matrix(A)
}
