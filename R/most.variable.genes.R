# Copyright (c) 2019 Leiden University Medical Center
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#     
#     The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#' Get top n most variable genes
#' 
#' Retrieve the list of genes identifiers for the genes which have the most 
#' variable expression in the given counts (based on voom normalized counts).
#' Rownames of the given count matrix are assumed to be gene identifiers.
#' @param counts A numeric matrix containing raw expression counts.
#' @param n The number of genes to be returned. Defaults to 500.
#'
#' @return A vector containing gene identifiers.
#' @export
#' @importFrom limma voom
#'
#' @examples
#' top500 <- most.variable.genes(counts)
most.variable.genes <- function(counts, n=500) {
    #TODO handle common datastructures besides data.frame: SE, RSE, MAE, etc.
    norm <- voom(as.matrix(counts), normalize.method="none")$E
    vars <- apply(norm, 1, var)
    top <- sort.list(vars, decreasing=TRUE)[1:n]
    rownames(counts[top,])
}

