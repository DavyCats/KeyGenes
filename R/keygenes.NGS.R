#' Title
#'
#' @param test blah
#' @param train blah
#' @param train.classes blah
#' @param genes blah
#' @param test.classes blah
#' @param verbose blah
#'
#' @return blah
#' @export
#'
#' @examples blah
keygenes.NGS <- function(test, train, train.classes, genes=NULL,
                        test.classes=NULL, verbose=FALSE) {
    if (extends(class(train), "SummarizedExperiment")) {
        suppressPackageStartupMessages(require(SummarizedExperiment))
        train.classes <- as.character(colData(train)[,train.classes])
        train <- assay(train)
    }
    
    if (extends(class(test), "SummarizedExperiment")) {
        suppressPackageStartupMessages(require(SummarizedExperiment))
        if (!is.null(test.classes)) test.classes <- as.character(
            colData(test)[,test.classes])
        test <- assay(test)
    }
    
    if (is.null(genes)) {
        if (verbose) message("Determining most variable genes")
        genes <- mostVariableGenes(train, n=500)
    }
    
    if (is.null(test.classes)) {
        if (verbose) message("Setting test classes to NA")
        test.classes <- rep(NA, times=ncol(test))
    }
    
    keygenes.NGS.run(as.matrix(test), as.matrix(train), train.classes, genes,
                     test.classes, verbose)
}


#' Title
#'
#' @param test blah
#' @param train blah
#' @param train.classes blah
#' @param genes blah
#' @param test.classes blah
#' @param verbose blah
#'
#' @return blah
#' @export
#' @importFrom limma voom
#' @importFrom glmnet cv.glmnet
#'
#' @examples blah
keygenes.NGS.run <- function(test, train, train.classes, genes,
                        test.classes, verbose=FALSE) {
    if (verbose) message("filtering and normalizing")
    # common genes
    common.genes <- intersect(row.names(train), row.names(test))
    
    # merge test and training sets and normalize together
    combined <- cbind(train[common.genes,], test[common.genes,])
    normalized <- voom(as.matrix(combined), normalize.method = "none")$E #TODO move normalization to keygenes.NGS? the rest here might be reusable for MA
    
    # filter for most variable genes
    norm.filtered <- t(normalized[rownames(normalized) %in% genes,])
    
    # split back up into train and test sets
    train.i <- 1:ncol(train) # indices for train samples
    test.i <- (1 + ncol(train)) : (ncol(train) + ncol(test)) # indices for test samples
    norm.train <- norm.filtered[train.i,]
    norm.test <- norm.filtered[test.i,]
    
    # duplicate samples
    norm.train.copy <- norm.train
    train.classes.copy <- train.classes
    while (any(table(train.classes) < 3)) {
        if (verbose) message("Not enough observations per class, duplicating data")
        norm.train <- rbind(data.frame(norm.train), 
                            data.frame(norm.train.copy))
        train.classes <- c(train.classes, train.classes.copy)
    }
    
    # get fold ids
    if (verbose) message("Determining folds")
    fold.id <- sample(rep(1:10, length(train.classes)), length(train.classes))
    while (any(sapply(1:10, 
                      function(i, f, m) {
                          any(table(m[f != i]) < 2)
                      },
                      fold.id, train.classes))) {
        if (verbose) message("...retrying")
        fold.id <- sample(rep(1:10, length(train.classes)), 
                          length(train.classes))
    }
    
    # create the model
    if (verbose) message("fitting model")
    cvfit <- cv.glmnet(as.matrix(norm.train), train.classes, 
                       family="multinomial", foldid = fold.id)
    
    # determine keyGenes(tm) for each class
    if (verbose) message("Retrieving key genes")
    coef <- coef(cvfit,  s=cvfit$lambda.min)
    class.genes <- lapply(coef, function(x){
        ind <- x@i[-1]+1
        x@Dimnames[[1]][ind]
    })
    
    # make prediction for the test set
    if (verbose) message("Making predicions")
    pred <- drop(predict(cvfit, newx=norm.test, type="response", 
                         s=cvfit$lambda.min))
    prediction.matrix <- t(pred)
    result <- data.frame(row.names=rownames(pred),
                         truth=test.classes,
                         predicted=colnames(pred)[apply(pred,1,which.max)])

    # Determine accuracy
    if (verbose) message("calculating accuracy")
    correct <- as.character(result$truth) == as.character(result$predicted)
    nCorrect <- sum(correct, na.rm = T)
    nTotal <- sum(!is.na(correct))
    accuracy <- nCorrect / nTotal
    
    # construct output
    if (verbose) message("preparing output")
    out <- new("KeyGenesResults",
               result=result,
               accuracy=accuracy,
               cvfit=cvfit,
               class.genes=class.genes,
               genes=genes,
               prediction.matrix=prediction.matrix,
               train=as.matrix(norm.train),
               train.classes=train.classes,
               test=as.matrix(norm.test),
               test.classes=test.classes)
    out
}
