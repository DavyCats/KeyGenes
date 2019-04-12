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
#' @slot duplicated.sample integer. The indexes of the samples which which were
#' duplicated in order to have a enough samples per class for training the model.
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

