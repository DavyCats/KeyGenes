context("keygenes NGS")

data("fetal_wo")
data("adult")

test_that("internal run function", {
    train <- SummarizedExperiment::assay(fetal_wo)
    test <- SummarizedExperiment::assay(adult)
    train_tissues <- as.character(SummarizedExperiment::colData(
        fetal_wo)$tissue)
    test_tissues <- as.character(SummarizedExperiment::colData(adult)$tissue)
    input_genes <- rownames(test)[1:500]
    
    result <- keygenes.NGS.run(test, train, train_tissues, input_genes,
                     test_tissues)
    # correct types
    expect(extends(class(result), "KeyGenesResults"))
    expect(extends(class(result@cvfit), "cv.glmnet"))
    
    # correct genes
    expect(all(sapply(result@class.genes, function(x){all(x %in% input_genes)})))
    expect(all(result@genes == input_genes))
    
    # make sure it doesn't resort to an intercept
    expect(all(apply(result@prediction.matrix, 1, function(x){min(x) != max(x)})))
    expect(all(apply(result@prediction.matrix, 2, function(x){min(x) != max(x)})))
    
    # correct tissues
    expect(all(result@train.classes == train_tissues))
    expect(all(result@test.classes == test_tissues))
    
    # TODO train and test should not equal input (should have been normalized/filtered)
})

test_that("SE input, no genes, no truth", {
    result <- keygenes.NGS(adult, fetal_wo, "tissue")
    
    # truth should be empty
    expect(all(is.na(result@result$truth)))
    
    # accuracy should be nan
    expect(is.nan(result@accuracy))
    
    # there should be 500 genes used.
    expect(length(result@genes) == 500)
}) 

test_that("SE input, no genes, with truth", {
    result <- keygenes.NGS(adult, fetal_wo, "tissue", test.classes = "tissue")
    expect(all(!is.na(result@result$truth)))
    expect(!is.nan(result@accuracy))
    
    #TODO test if tissue/classes are correct
}) 

test_that("SE input, with genes, with truth", {
    result <- keygenes.NGS(adult, fetal_wo, "tissue", 
                         genes = rownames(fetal_wo)[1:500], 
                         test.classes = "tissue")
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

