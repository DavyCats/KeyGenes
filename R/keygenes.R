setClassUnion("OptionalChar",   c("character",   "logical"))

#' Title
#'
#' @slot result data.frame. 
#' @slot accuracy numeric. 
#' @slot cvfit ANY. 
#' @slot class.genes list. 
#' @slot prediction.matrix matrix. 
#' @slot train matrix. 
#' @slot train.classes character. 
#' @slot test matrix. 
#' @slot test.classes OptionalChar. 
#' @slot genes character. 
#'
#' @return blah
#' @export
#' @import glmnet
#'
#' @examples blah
setClass("KeyGenesResults", slots=list(result="data.frame",
                                       accuracy="numeric",
                                       cvfit="ANY", #tested for correct type below
                                       class.genes="list",
                                       prediction.matrix="matrix",
                                       train="matrix",
                                       train.classes="character",
                                       test="matrix",
                                       test.classes="OptionalChar",
                                       genes="character"),
         validity = function(object){
             extends(class(object@cvfit), "cv.glmnet")
         })

