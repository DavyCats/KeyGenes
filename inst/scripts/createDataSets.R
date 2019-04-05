hESC <- "sourceData/hESC/test_forster.txt"
hESC.data <- read.table(hESC, sep="\t", header=TRUE, row.names=1,
                        check.names=FALSE)[,c("DH3.hES_Stem Cells",
                                              "DH6.hES_Stem Cells",
                                              "DH7.hES_Stem Cells")]
colnames(hESC.data) <- c("hESC Stem Cells.DH3", "hESC Stem Cells.DH6",
                         "hESC Stem Cells.DH7")
fetal <- "sourceData/human_fetal/training_fetal.txt"
fetal.data <- read.table(fetal, sep="\t", header=TRUE, row.names=1,
                         check.names=FALSE)
adult <- "sourceData/human_adult/training_adult.txt"
adult.data <- read.table(adult, sep="\t", header=TRUE, row.names=1,
                         check.names=FALSE)
colnames(adult.data) <- make.names(colnames(adult.data), unique = T)

all.genes <- union(rownames(hESC.data), rownames(fetal.data))
all.genes <- union(all.genes, rownames(adult.data))

missing.hESC <- setdiff(all.genes, row.names(hESC.data))
missing.fetal <- setdiff(all.genes, rownames(fetal.data))
missing.adult <- setdiff(all.genes, row.names(adult.data))

hESC.data[missing.hESC,] <- NA
fetal.data[missing.fetal,] <- NA
adult.data[missing.adult,] <- NA

setdiff(all.genes, row.names(hESC.data))
setdiff(all.genes, rownames(fetal.data))
setdiff(all.genes, row.names(adult.data))

all.counts <- cbind(hESC.data, fetal.data)
all.counts <- cbind(all.counts, adult.data)

coldata <- data.frame(row.names = colnames(all.counts))

# age
coldata[grepl("_9", colnames(all.counts)), "age"] <- "1st trimester"
coldata[grepl("_(22|16-18)", colnames(all.counts)), "age"] <- "2nd trimester"
coldata[grepl("_adult", colnames(all.counts)), "age"] <- "Adult"
coldata[grepl("hESC Stem Cells", colnames(all.counts)), "age"] <- "hESC stem cells"

# organ
organs <- sapply(strsplit(colnames(all.counts), "_"), "[[", 1)
organs <- sapply(strsplit(organs, "\\."), "[[", 1)
coldata[, "organ"] <- organs

training.data <- SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=as.matrix(all.counts)),
    colData=coldata)

save(training.data, file="data/training_data.RData", compress = "xz", compression_level = 9)

