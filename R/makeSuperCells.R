#' Create super-cells by aggregating similar cells
#' @description Constructs metacells/super-cells from single-cell data using 
#'   graph-based clustering in PCA space
#' @param X Gene expression matrix (genes x cells)
#' @param genes.use Genes to use for PCA
#' @param genes.exclude Genes to exclude from analysis
#' @param n.var.genes Number of variable genes to select
#' @param gamma Graining level (ratio of cells to super-cells)
#' @param k.knn Number of nearest neighbors for graph construction
#' @param do.scale Whether to scale data before PCA
#' @param n.pc Number of principal components
#' @param fast.pca Use truncated SVD for speed
#' @param do.approx Use approximate algorithm for large datasets
#' @param approx.N Sample size for approximation
#' @param directed Build directed kNN graph
#' @param use.nn2 Use RANN::nn2 for kNN
#' @param seed Random seed
#' @param igraph.clustering Clustering algorithm
#' @param return.singlecell.NW Return single-cell network
#' @param return.hierarchical.structure Return hierarchy
#' @param block.size Block size for large data processing
#' @param weights Cell weights for aggregation
#' @param do.median.norm Normalize by median
#' @return Aggregated gene expression matrix (genes x super-cells)
#' @importFrom stats aggregate as.dist cor median prcomp sd var
#' @importFrom irlba irlba
#' @importFrom igraph cluster_walktrap cut_at cluster_louvain contract simplify
#' @import Matrix
#' @noRd
makeSuperCells <- function(X,
                           genes.use = NULL,
                           genes.exclude = NULL,
                           n.var.genes = min(1000, nrow(X)),
                           gamma = 10,
                           k.knn = 5,
                           do.scale = TRUE,
                           n.pc = 25,
                           fast.pca = TRUE,
                           do.approx = FALSE,
                           approx.N = 20000,
                           directed = FALSE,
                           use.nn2 = TRUE,
                           seed = 12345,
                           igraph.clustering = c("walktrap", "louvain"),
                           return.singlecell.NW = TRUE,
                           return.hierarchical.structure = TRUE,
                           block.size = 10000,
                           weights = NULL,
                           do.median.norm = FALSE) {
  N.c <- ncol(X)
  GE <- X

  if (is.null(rownames(X))) {
    if (!is.null(genes.use) || !is.null(genes.exclude)) {
      stop(
        "rownames(X) is NULL but genes.use or genes.exclude is specified.\n",
        "Gene expression matrix X must have genes as rownames."
      )
    } else {
      warning(
        "rownames(X) is NULL.\n",
        "Gene expression matrix X is expected to have genes as rownames.\n",
        "Gene names will be created automatically as 'gene_1', 'gene_2', etc."
      )
      rownames(X) <- paste("gene", seq_len(nrow(X)), sep = "_")
    }
  }

  if (is.null(colnames(X))) {
    warning(
      "colnames(X) is NULL.\n",
      "Gene expression matrix X is expected to have cell IDs as colnames.\n",
      "Cell IDs will be created automatically as 'cell_1', 'cell_2', etc."
    )
    colnames(X) <- paste("cell", seq_len(N.c), sep = "_")
  }

  X <- X[rowSums(X) > 0, ]
  keep.genes <- setdiff(rownames(X), genes.exclude)
  X <- X[keep.genes, ]


  if (is.null(genes.use)) {
    n.var.genes <- min(n.var.genes, nrow(X))
    if (N.c > 50000) {
      set.seed(seed)
      idx <- sample(N.c, 50000)
      X_sub <- X[, idx]
      row_means <- rowMeans(X_sub)
      gene.var <- rowMeans((X_sub - row_means)^2) * ncol(X_sub) / (ncol(X_sub) - 1)
    } else {
      row_means <- rowMeans(X)
      gene.var <- rowMeans((X - row_means)^2) * N.c / (N.c - 1)
    }
    names(gene.var) <- rownames(X)
    genes.use <- names(sort(gene.var, decreasing = TRUE))[1:n.var.genes]
  }

  if (length(intersect(genes.use, genes.exclude)) > 0) {
    stop("Sets of genes.use and genes.exclude have non-empty intersection")
  }

  genes.use <- genes.use[genes.use %in% rownames(X)]
  X <- X[genes.use, ]

  if (do.approx & approx.N >= N.c) {
    do.approx <- FALSE
    warning(
      "N.approx is larger or equal to the number of single cells, thus, an exact simplification will be performed"
    )
  }

  if (do.approx & (approx.N < round(N.c / gamma))) {
    approx.N <- round(N.c / gamma)
    warning(paste("N.approx is set to N.SC", approx.N))
  }

  if (do.approx & ((N.c / gamma) > (approx.N / 3))) {
    warning(
      "N.approx is not much larger than desired number of super-cells, so an approximate simplification may take longer than an exact one!"
    )
  }

  if (do.approx) {
    set.seed(seed)
    approx.N <- min(approx.N, N.c)
    presample <-
      sample(1:N.c, size = approx.N, replace = FALSE)
    presampled.cell.ids <- colnames(X)[sort(presample)]
    rest.cell.ids <- setdiff(colnames(X), presampled.cell.ids)
  } else {
    presampled.cell.ids <- colnames(X)
    rest.cell.ids <- c()
  }

  X.for.pca <-
    Matrix::t(X[genes.use, presampled.cell.ids])
  if (do.scale) {
    X.for.pca <- scale(X.for.pca)
  }
  X.for.pca[is.na(X.for.pca)] <- 0

  if (is.null(n.pc[1]) |
    min(n.pc) < 1) {
    stop("Please, provide a range or a number of components to use: n.pc")
  }
  if (length(n.pc) == 1) {
    n.pc <- 1:n.pc
  }

  if (fast.pca & (N.c < 1000)) {
    # warning("Normal PCA is computed because number of cell is low for irlba::irlba()")
    fast.pca <- FALSE
  }

  if (!fast.pca) {
    PCA.presampled <-
      prcomp(
        X.for.pca,
        rank. = max(n.pc),
        scale. = F,
        center = F
      )
  } else {
    PCA.presampled <-
      irlba::irlba(X.for.pca, max(n.pc))
    PCA.presampled$x <-
      PCA.presampled$u %*% diag(PCA.presampled$d)
    PCA.presampled$rotation <- PCA.presampled$v
  }

  sc.nw <-
    buildKNN(
      X = PCA.presampled$x[, n.pc],
      k = k.knn,
      from = "coordinates",
      use.nn2 = use.nn2,
      dist_method = "euclidean",
      directed = directed
    )

  # simplify

  k <- round(N.c / gamma)

  if (igraph.clustering[1] == "walktrap") {
    g.s <- igraph::cluster_walktrap(sc.nw$graph.knn)
    g.s$membership <- igraph::cut_at(g.s, k)
  } else if (igraph.clustering[1] == "louvain") {
    warning(paste(
      "igraph.clustering =",
      igraph.clustering,
      ", gamma is ignored"
    ))
    g.s <- igraph::cluster_louvain(sc.nw$graph.knn)
  } else {
    stop(
      paste(
        "Unknown clustering method (",
        igraph.clustering,
        "), please use louvain or walktrap"
      )
    )
  }

  membership.presampled <- g.s$membership
  names(membership.presampled) <- presampled.cell.ids

  SC.NW <-
    igraph::contract(sc.nw$graph.knn, membership.presampled)
  SC.NW <-
    igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb = "sum")

  if (do.approx) {
    # Approximate mode not fully implemented yet
    # Using presampled membership as fallback
    membership.all <- membership.presampled[colnames(X)]
  } else {
    membership.all <- membership.presampled[colnames(X)]
  }

  membership <- membership.all
  X <- GE

  N.SC <- max(membership)
  supercell_size <- base::as.vector(table(membership))
  j <-
    rep(1:N.SC, supercell_size) # column indices of matrix M.AV that, whene GE.SC <- ge %M.AV%

  goups.idx <- base::split(seq_len(ncol(X)), membership) # plyr::split_indices(membership)
  i <-
    unlist(goups.idx) # row indices of matrix M.AV that, whene GE.SC <- ge %M.AV%

  if (is.null(weights)) {
    M.AV <- Matrix::sparseMatrix(i = i, j = j)
    GE <- X %*% M.AV
    GE <- sweep(GE, 2, supercell_size, "/")
  } else {
    if (length(weights) != length(membership)) {
      stop("weights must be the same length as groups or NULL in case of unweighted averaging")
    }
    M.AV <- Matrix::sparseMatrix(i = i, j = j, x = weights[i])
    GE <- GE %*% M.AV

    weighted_supercell_size <-
      unlist(lapply(
        goups.idx,
        FUN = function(x) {
          sum(weights[x])
        }
      ))
    GE <- sweep(GE, 2, weighted_supercell_size, "/")
  }

  if (do.median.norm) {
    GE_shifted <- GE + 0.01
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      row_medians <- matrixStats::rowMedians(as.matrix(GE_shifted))
    } else {
      row_medians <- apply(as.matrix(GE_shifted), 1, median)
    }
    GE <- GE_shifted / row_medians
  }

  return(GE)
}
