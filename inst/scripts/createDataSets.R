hESC <- "sourceData/hESC/test_forster.txt"
hESC.data <- read.table(hESC, sep="\t", header=TRUE, row.names=1,
                        check.names=FALSE)[,c("DH3.hES_Stem Cells",
                                              "DH6.hES_Stem Cells",
                                              "DH7.hES_Stem Cells")]
colnames(hESC.data) <- c("hESC Stem Cells.DH3", "hESC Stem Cells.DH6",
                         "hESC Stem Cells.DH7")
hESC2 <- "sourceData/hESC/2013-11-04_combined_all_data.csv"
hESC2.data <- read.table(hESC2, sep=",", header=TRUE, row.names=1,
                         check.names=FALSE)[, "set_3_B6_GCCAAT_L005", drop=F]
colnames(hESC2.data) <- "hESC Stem Cells.1"
fetal <- "sourceData/human_fetal/training_fetal.txt"
fetal.data <- read.table(fetal, sep="\t", header=TRUE, row.names=1,
                         check.names=FALSE)
adult <- "sourceData/human_adult/training_adult.txt"
adult.data <- read.table(adult, sep="\t", header=TRUE, row.names=1,
                         check.names=FALSE)
colnames(adult.data) <- make.names(colnames(adult.data), unique = T)

all.genes <- union(rownames(hESC.data), rownames(hESC2.data))
all.genes <- union(all.genes, rownames(fetal.data))
all.genes <- union(all.genes, rownames(adult.data))

missing.hESC <- setdiff(all.genes, row.names(hESC.data))
missing.hESC2 <- setdiff(all.genes, row.names(hESC2.data))
missing.fetal <- setdiff(all.genes, rownames(fetal.data))
missing.adult <- setdiff(all.genes, row.names(adult.data))

hESC.data[missing.hESC,] <- NA
hESC2.data[missing.hESC2,] <- NA
fetal.data[missing.fetal,] <- NA
adult.data[missing.adult,] <- NA

setdiff(all.genes, row.names(hESC.data))
setdiff(all.genes, row.names(hESC2.data))
setdiff(all.genes, rownames(fetal.data))
setdiff(all.genes, row.names(adult.data))

all.counts <- cbind(hESC.data, hESC2.data)
all.counts <- cbind(all.counts, fetal.data)
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

save(training.data, file="data/training.data.RData", compress = "xz", compression_level = 9)

