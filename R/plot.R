#' Title
#'
#' @param data blah
#' @param clusterTissues blah
#' @param clusterSamples  blah
#'
#' @return blah
#' @export
#'
#' @examples blah
keygenes.heatmap <- function(data, clusterTissues=F, clusterSamples=F){
    require(ggplot2)
    
    if (!extends(class(data), "KeyGenesResults")){
        stop("data needs to be a KeyGenesResults object!")
    }
    
    pred.final <- data@prediction.matrix
    
    if (clusterTissues || clusterSamples) scaled <- scale(data@prediction.matrix)
    if (clusterTissues) {
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
