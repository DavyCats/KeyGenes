#' Get top n most variable genes
#' 
#' Retrieve the list of genes identifiers for the genes which have the most 
#' variable expression in the given counts. Rownames of the given count matrix
#' are assumed to be gene identifiers.
#' @param counts A numeric matrix containing raw expression counts.
#' @param n The number of genes to be returned. Defaults to 500.
#'
#' @return A vector containing gene identifiers.
#' @export
#' @importFrom limma voom
#'
#' @examples
#' top500 <- mostVariableGenes(counts)
mostVariableGenes <- function(counts, n=500) {
    #TODO handle common datastructures besides data.frame: SE, RSE, MAE, etc.
    norm <- voom(as.matrix(counts), normalize.method="none")$E
    vars <- apply(norm, 1, var)
    top <- sort.list(vars, decreasing=TRUE)[1:n]
    rownames(counts[top,])
}

