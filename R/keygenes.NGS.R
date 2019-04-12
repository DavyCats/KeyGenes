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

#' Predict the class (eg. organ) for NGS (RNA-seq) samples
#'
#' Predict the class (eg. organ or maturation) for the test samples based on
#' a set of training samples, using RNA-seq expression measures.
#'
#' @param test The counts table of the test data or a SummarizedExperiment 
#' containing it.
#' @param train The counts table of the training data or a SummarizedExperiment
#' containing it.
#' @param train.classes The class each sample in the training set belongs to.
#' Each sample from the test set will be assigned to one of these. Normally a
#' character vector with one item per train sample, if `train` is a 
#' SummarizedExperiment the colData column containing the classes.
#' @param genes The genes to be used for training and testing. If NULL
#' \link[KeyGenes]{most.variable.genes} will be called wih n = 500. Defaults to
#' NULL
#' @param test.classes See train.classes, but for the test samples. Used to
#' calculate the accuracy of the predictions. If NULL the accuracy will be NaN.
#' Defaults to NULL 
#' @param verbose Should progress be reported, defaults to FALSE
#'
#' @return A \link[KeyGenes]{KeyGenesResults} object containing the prediction 
#' results.
#' @export
#'
#' @examples 
#' data(adult)
#' data(fetal_wo)
#' result <- keygenes.NGS(adult, fetal_wo, "tissue")
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
    
    if (verbose) message("Selecting shared genes")
    test <- as.matrix(test)
    train <- as.matrix(train)
    common.genes <- intersect(row.names(train), row.names(test))
    test <- test[common.genes,]
    train <- train[common.genes,]
    
    # drop classes with only one sample
    drop <- names(which(table(train.classes) < 2))
    for (x in drop) warning("Not enough samples in training data for class '",
                            x, "', will be dropped")
    train <- train[,! train.classes %in% drop]
    train.classes <- train.classes[! train.classes %in% drop]
    
    if (is.null(genes)) {
        if (verbose) message("Determining most variable genes")
        genes <- most.variable.genes(train, n=500)
    }
    
    if (is.null(test.classes)) {
        if (verbose) message("Setting test classes to NA")
        test.classes <- rep(NA, times=ncol(test))
    }
    
    keygenes.NGS.run(test, train, train.classes, genes,
                     test.classes, verbose)
}


#' Predict the class (eg. organ) for NGS (RNA-seq) samples
#'
#' Predict the class (eg. organ or maturation) for the test samples based on
#' a set of training samples, using RNA-seq expression measures.
#'
#' @param test A matrix of expression measures for the test set
#' @param train A matrix of expression measures for the training set
#' @param train.classes A vector with the classes of the training samples. Each
#' sample from the test set will be assigned to one of these
#' @param genes The genes to be used for training and testing
#' @param test.classes A vector with the classes of the test samples. Used to
#' calculate accuracy. May be NULL
#' @param verbose Should progress be reported, defaults to FALSE
#'
#' @return A  \link[KeyGenes]{KeyGenesResults} object containing the prediction 
#' results.
#' @importFrom limma voom
#' @importFrom glmnet cv.glmnet
keygenes.NGS.run <- function(test, train, train.classes, genes,
                        test.classes, verbose=FALSE) {
    if (verbose) message("Normalizing and filtering")
    # merge test and training sets and normalize together
    combined <- cbind(train, test)
    normalized <- voom(as.matrix(combined), normalize.method = "none")$E
    
    # filter for most variable genes
    norm.filtered <- t(normalized[rownames(normalized) %in% genes,])
    
    # split back up into train and test sets
    train.i <- 1:ncol(train) # indices for train samples
    test.i <- (1 + ncol(train)) : (ncol(train) + ncol(test)) # indices for test samples
    norm.train <- norm.filtered[train.i,]
    norm.test <- norm.filtered[test.i,]
    
    # determine folds
    if (verbose) message("Determining folds")
    n.folds <- min(c(10, length(train.classes)))
    if (verbose) message("Number of folds: ", n.folds)
    lacking.classes <- names(which(table(train.classes) < 3))
    not.lacking.classes <- names(which(table(train.classes) > 2))
    fold.id <- sample(rep(1:n.folds, length(train.classes)), 
                      length(train.classes))
    while (any(sapply(1:n.folds, 
                      function(i, f, m, l, n) {
                          not.in.fold <- table(m[f != i])
                          if ( ! all(train.classes %in% names(not.in.fold))) 
                              return(T)
                          if ( length(m[f==i]) == 0 )
                              return(T)
                          any(not.in.fold[n] < 2) | any(not.in.fold[l] < 1)
                      },
                      fold.id, train.classes, lacking.classes,
                      not.lacking.classes))) {
        if (verbose) message("...problematic class distribution, retrying")
        fold.id <- sample(rep(1:n.folds, length(train.classes)), 
                          length(train.classes))
    }
    
    # duplicate classes with only two samples
    if (verbose && length(lacking.classes) > 0)
        message("Not enough observations per class, duplicating data")
    
    duplication.samples <- which(train.classes %in% lacking.classes)
    duplication <- norm.train[duplication.samples,]
    duplication.classes <- train.classes[duplication.samples]
    duplication.fold.id <- fold.id[duplication.samples]
    
    final.train <- rbind(norm.train, duplication)
    final.train.classes <- c(train.classes , duplication.classes)
    final.fold.id <- c(fold.id, duplication.fold.id)
    train.weights <- c(rep(1, length(train.classes)), rep(0.5, length(duplication.classes)))
    train.weights[duplication.samples] <- 0.5
    
    # create the model
    if (verbose) message("fitting model")
    cvfit <- cv.glmnet(as.matrix(final.train), final.train.classes,
                       nfolds=n.folds, family="multinomial",
                       foldid = final.fold.id)
    
    # determine keyGenes for each class
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
               test.classes=test.classes,
               duplicated.samples=duplication.samples)
    out
}
