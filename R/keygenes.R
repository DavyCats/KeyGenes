setClassUnion("OptionalChar",   c("character",   "logical"))

#' KeyGenes prediction results
#'
#' @slot result data.frame. The predictions.
#' @slot accuracy numeric.  Accuracy of the predictions or NaN.
#' @slot cvfit ANY. The model used for the predictions.
#' @slot class.genes list. The key genes for each class.
#' @slot prediction.matrix matrix. The raw prediction values.
#' @slot train matrix. The final effective training set.
#' @slot train.classes character. The classes for each sample in the training set.
#' @slot test matrix. The final effective test set.
#' @slot test.classes OptionalChar.  The classes for the test set.
#' @slot genes character. The genes used for training and prediction.
#'
#' @export
#' @import glmnet
setClass("KeyGenesResults", slots=list(result="data.frame",
                                       accuracy="numeric",
                                       cvfit="ANY", #tested for correct type below
                                       class.genes="list",
                                       prediction.matrix="matrix",
                                       train="matrix",
                                       train.classes="character",
                                       test="matrix",
                                       test.classes="OptionalChar",
                                       genes="character",
                                       duplicated.samples="integer"),
         validity = function(object){
             extends(class(object@cvfit), "cv.glmnet")
         })

