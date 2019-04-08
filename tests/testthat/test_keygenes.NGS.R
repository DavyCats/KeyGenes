context("keygenes NGS")

data("training.data")
training.data <- training.data[complete.cases(SummarizedExperiment::assay(training.data)),]
train.se <- training.data[,training.data$age == "Adult"]
train.se <- train.se[, ! train.se$organ %in% c("muscle", "skin")]
test.se <- training.data[,training.data$age %in% c("1st trimester", "2nd trimester")]
train <- SummarizedExperiment::assay(train.se)
test <- SummarizedExperiment::assay(test.se)

test_that("internal run function", {
    train.organs <- as.character(SummarizedExperiment::colData(train.se)$organ)
    test.organs <- as.character(SummarizedExperiment::colData(test.se)$organ)
    
    common.genes <- intersect(rownames(train), row.names(test))
    
    result <- keygenes.NGS.run(test[common.genes,1:5],
                               train[common.genes,], 
                               train.organs, common.genes[1:500], test.organs[1:5])
    
    # correct types
    expect(extends(class(result), "KeyGenesResults"))
    expect(extends(class(result@cvfit), "cv.glmnet"))
    
    # correct genes
    expect(all(sapply(result@class.genes, function(x){all(x %in% common.genes[1:500])})))
    expect(all(result@genes == common.genes[1:500]))
    
    # make sure it doesn't resort to an intercept
    expect(all(apply(result@prediction.matrix, 1, function(x){min(x) != max(x)})))
    expect(all(apply(result@prediction.matrix, 2, function(x){min(x) != max(x)})))
    
    # correct organs
    expect(all(result@train.classes == train.organs))
    expect(all(result@test.classes == test.organs[1:5]))
    
    # TODO train and test should not equal input (should have been normalized/filtered)
    # TODO duplicated samples
})

test_that("SE input, no genes, no truth", {
    result <- keygenes.NGS(test.se, train.se, "organ")
    
    # truth should be empty
    expect(all(is.na(result@result$truth)))
    
    # accuracy should be nan
    expect(is.nan(result@accuracy))
    
    # there should be 500 genes used.
    expect(length(result@genes) == 500)
}) 

test_that("SE input, no genes, with truth", {
    result <- keygenes.NGS(test.se, train.se, "organ", test.classes = "organ")
    expect(all(!is.na(result@result$truth)))
    expect(!is.nan(result@accuracy))
    
    #TODO test if tissue/classes are correct
}) 

test_that("SE input, with genes, with truth", {
    result <- keygenes.NGS(test.se, train.se, "organ", 
                         genes = rownames(train.se)[1:500], 
                         test.classes = "organ")
    expect(all(!is.na(result@result$truth)))
    expect(!is.nan(result@accuracy))
    
    #TODO check if genes are correct
}) 

#TODO
# test non SE input
    # with truth 
        # check that classes (train.classes and test.classes) are correct
    # without truth
    # genes input already tested no need to do so again

