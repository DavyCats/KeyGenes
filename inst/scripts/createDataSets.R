T1 <- "sourceData/human_fetal_wo_1t/training_fetal_wo_1t.txt"
T1Data <- read.table(T1, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
T2 <- "sourceData/human_fetal_wo_2t/training_fetal_wo_2t.txt"
T2Data <- read.table(T2, sep="\t", header=TRUE, row.names=1, check.names=FALSE)

samplesT1 <- colnames(T1Data)
tissuesT1 <- sapply(strsplit(samplesT1, "_"), "[", 1)
samplesT2 <- colnames(T2Data)
tissuesT2 <- sapply(strsplit(samplesT2, "_"), "[", 1)

fetalCol <- data.frame(row.names = c(colnames(T1Data), colnames(T2Data)))
fetalCol$tissue <- c(tissuesT1, tissuesT2)
fetalCol$trimester <- c(rep(1, times=length(tissuesT1)), 
                        rep(2, times=length(tissuesT2)))

trainData <- cbind(T1Data, T2Data)

fetal_wo <- SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=as.matrix(trainData)),
    colData=fetalCol)

save(fetal_wo, file="data/fetal_wo.RData", compress = "xz", compression_level = 9)

# ---

adultPath <- "sourceData/human_adult/training_adult.txt"
adultData <- read.table(adult, sep="\t", header=TRUE, row.names=1, 
                        check.names=FALSE)
colnames(adultData) <- make.names(colnames(adultData), unique = T)

adultTissues <- sapply(strsplit(colnames(adultData), "_"), "[", 1)
adultCol <- data.frame(row.names = c(colnames(adultData)), tissue=adultTissues)

adult <- SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=as.matrix(adultData)),
    colData=adultCol)

save(adult, file="data/adult.RData", compress = "xz", compression_level = 9)
