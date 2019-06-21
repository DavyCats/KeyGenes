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


hESC <- "sourceData/hESC/test_forster.txt"
hESC.data <- read.table(hESC, sep="\t", header=TRUE, row.names=1,
                        check.names=FALSE)[,c("DH3.hES_Stem Cells",
                                              "DH6.hES_Stem Cells",
                                              "DH7.hES_Stem Cells")]
colnames(hESC.data) <- c("hESC stem cells.DH3", "hESC stem cells.DH6",
                         "hESC stem cells.DH7")
hESC2 <- "sourceData/hESC/2013-11-04_combined_all_data.csv"
hESC2.data <- read.table(hESC2, sep=",", header=TRUE, row.names=1,
                         check.names=FALSE)[, "set_3_B6_GCCAAT_L005", drop=F]
colnames(hESC2.data) <- "hESC stem cells.1"
hESC2.data <- hESC2.data[! grepl("LRG_[0-9]+", row.names(hESC2.data)),, drop=F]
fetal <- "sourceData/human_fetal/training_fetal.txt"
fetal.data <- read.table(fetal, sep="\t", header=TRUE, row.names=1,
                         check.names=FALSE)
adult <- "sourceData/human_adult/training_adult.txt"
adult.data <- read.table(adult, sep="\t", header=TRUE, row.names=1,
                         check.names=FALSE)
colnames(adult.data) <- make.names(colnames(adult.data), unique = T)

all.genes <- union(row.names(hESC.data), row.names(hESC2.data))
all.genes <- union(all.genes, row.names(fetal.data))
all.genes <- union(all.genes, row.names(adult.data))

missing.hESC <- setdiff(all.genes, row.names(hESC.data))
missing.hESC2 <- setdiff(all.genes, row.names(hESC2.data))
missing.fetal <- setdiff(all.genes, row.names(fetal.data))
missing.adult <- setdiff(all.genes, row.names(adult.data))

hESC.data[missing.hESC,] <- NA
hESC2.data[missing.hESC2,] <- NA
fetal.data[missing.fetal,] <- NA
adult.data[missing.adult,] <- NA

setdiff(all.genes, row.names(hESC.data))
setdiff(all.genes, row.names(hESC2.data))
setdiff(all.genes, row.names(fetal.data))
setdiff(all.genes, row.names(adult.data))

all.counts <- merge(hESC.data, hESC2.data, by=0)
all.counts <- merge(all.counts, fetal.data, by.x=1, by.y=0)
all.counts <- merge(all.counts, adult.data, by.x=1, by.y=0)
row.names(all.counts) <- all.counts$Row.names
all.counts <- all.counts[,colnames(all.counts) != "Row.names"]

coldata <- data.frame(row.names = colnames(all.counts))

# age
coldata[grepl("_9", colnames(all.counts)), "age"] <- "1st trimester"
coldata[grepl("_(22|16-18)", colnames(all.counts)), "age"] <- "2nd trimester"
coldata[grepl("_adult", colnames(all.counts)), "age"] <- "Adult"
coldata[grepl("hESC stem cells", colnames(all.counts)), "age"] <- "hESC stem cells"

# organ
organs <- sapply(strsplit(colnames(all.counts), "_"), "[[", 1)
organs <- sapply(strsplit(organs, "\\."), "[[", 1)
coldata[, "organ"] <- organs

rowdata <- data.frame(row.names = row.names(all.counts))
ensembl <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensembl_ids <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_id", "hgnc_symbol"), mart=ensembl)

# hgnc_id
hgnc_id <- aggregate(hgnc_id ~ ensembl_gene_id, data=ensembl_ids, FUN=paste, collapse=";")
row.names(hgnc_id) <- hgnc_id$ensembl_gene_id
rowdata[,"hgnc_id"] <- hgnc_id[row.names(rowdata),"hgnc_id"]

# symbol
symbol <- aggregate(hgnc_symbol ~ ensembl_gene_id, data=ensembl_ids, FUN=paste, collapse=";")
row.names(symbol) <- symbol$ensembl_gene_id
rowdata[,"symbol"] <- symbol[row.names(rowdata),"hgnc_symbol"]

rowdata[rowdata == ""] <- NA

training.data <- SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=as.matrix(all.counts)),
    colData=coldata, rowData=rowdata)

save(training.data, file="data/training.data.RData", compress = "xz", compression_level = 9)

