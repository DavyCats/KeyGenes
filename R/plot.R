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

#' Keygenes heatmap
#' 
#' Draw a heatmap for KeyGenes results.
#'
#' @param data A KeyGenesResults containing prediction results.
#' @param clusterClasses Whether or not to cluster the classes. Defaults to FALSE.
#' @param clusterSamples Whether or not to cluster the samples. Defaults to FALSE.
#'
#' @return A ggplot displaying a heatmap of the predictions.
#' @export
#'
#' @examples 
#' data(adult)
#' data(fetal_wo)
#' result <- keygenes.NGS(adult, fetal_wo, "tissue")
#' keygenes.heatmap(result)
keygenes.heatmap <- function(data, clusterClasses=F, clusterSamples=F){
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(reshape2))
    
    if (!extends(class(data), "KeyGenesResults")){
        stop("data needs to be a KeyGenesResults object!")
    }
    
    pred.final <- data@prediction.matrix
    
    if (clusterClasses || clusterSamples) scaled <- scale(data@prediction.matrix)
    if (clusterClasses) {
        ord <- hclust( dist(scaled, method = "euclidean"), method = "ward.D" )$order
        pred.final <- pred.final[ord,]
    }
    if (clusterSamples) {
        ord <- hclust( dist(t(scaled), method = "euclidean"), method = "ward.D" )$order # samples
        pred.final <- pred.final[,ord]
    }
    
    meltedPrediction <- melt(pred.final, varnames = c("tissue", "sample"))
    g <- ggplot(data=meltedPrediction)
    g <- g + theme_minimal()
    g <- g + geom_tile(aes(fill=value, x=sample, y=tissue), color="gray")
    g <- g + scale_fill_gradient2(low="black", mid="white", high = "green2", 
                                  midpoint = 0.5)
    g <- g + labs(x="", y="")
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                   axis.ticks = element_blank())
    g
}

#' KeyGenes dendrogram
#' 
#' Draw a dendrogram for KeyGenes results.
#'
#' @param data A KeyGenesResults containing prediction results.
#'
#' @return A ggplot displaying a dendrogram of the training and test data.
#' @export
#'
#' @examples
#' #' data(adult)
#' data(fetal_wo)
#' result <- keygenes.NGS(adult, fetal_wo, "tissue")
#' keygenes.dendrogram(result)
keygenes.dendrogram <- function(data) {
    suppressPackageStartupMessages(library(ggdendro))
    
    classes <- data@train.classes
    uniqueClasses <- unique(classes)
    
    genes <- unique(unlist(data@class.genes))
    
    meansPerClass <- apply(data@train, 2, function(x) {
        out <- c(rep(0, times=length(uniqueClasses)))
        names(out) <- uniqueClasses
        for (y in uniqueClasses) {
            out[y] <- mean(x[classes == y])
        }
        out
    })
    meansPerClass <- t(meansPerClass)
    
    countsAndMeans <- cbind(meansPerClass[genes,], t(data@test)[genes,])
    #countsAndMeans <- voom(as.matrix(countsAndMeans))$E
    
    corr <- cor(countsAndMeans, method="pearson")
    dis <- dist(corr)
    hc <- hclust(dis)
    ggdendrogram(hc, rotate=T)
}
