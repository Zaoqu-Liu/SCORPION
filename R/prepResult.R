#' Prepare PANDA output results
#' @param zScale Whether to use Z-scaling
#' @param output Output networks to include
#' @param regulatoryNetwork The regulatory network matrix
#' @param geneCoreg The gene co-regulation network matrix
#' @param tfCoopNetwork The TF cooperation network matrix
#' @param edgelist Whether to return as edgelist
#' @param motif The motif data
#' @return A list containing network matrices and statistics
#' @noRd
prepResult <- function(zScale, output, regulatoryNetwork, geneCoreg, tfCoopNetwork, edgelist, motif) {
  numGenes <- dim(geneCoreg)[1]
  numTFs <- dim(tfCoopNetwork)[1]
  numEdges <- sum(regulatoryNetwork != 0)
  list(
    regNet = regulatoryNetwork,
    coregNet = geneCoreg,
    coopNet = tfCoopNetwork,
    numGenes = numGenes,
    numTFs = numTFs,
    numEdges = numEdges
  )
}
