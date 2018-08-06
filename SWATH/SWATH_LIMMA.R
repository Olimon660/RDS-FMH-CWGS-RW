library(limma)
library(qvalue)
library(Biobase)
source("SWATH/SWATH_LIMMA_TOOLS.R")

noRep <- function(x) { gsub('_r1', '', gsub('_r2', '', gsub('_r3', '', x))) }

data <- read.table("SWATH/data2.tsv", row.names=1, header=TRUE, sep='\t')

stopifnot(is.na(data[row.names(data) == "sp|Q9C0F1|CEP44_HUMAN",]$WI38_r1))
stopifnot(is.na(data[row.names(data) == "sp|P40305|IFI27_HUMAN",]$IIICF_P7_r2))
stopifnot(is.na(data[row.names(data) == "sp|P08047|SP1_HUMAN",]$LFS_05F_24_pre_r3))

data[row.names(data) == "sp|Q9C0F1|CEP44_HUMAN",]$WI38_r1 <- mean(c(data[row.names(data) == "sp|Q9C0F1|CEP44_HUMAN",]$WI38_r2, data[row.names(data) == "sp|Q9C0F1|CEP44_HUMAN",]$WI38_r3))
data[row.names(data) == "sp|P40305|IFI27_HUMAN",]$IIICF_P7_r2 <- mean(c(data[row.names(data) == "sp|P40305|IFI27_HUMAN",]$IIICF_P7_r1, data[row.names(data) == "sp|P40305|IFI27_HUMAN",]$IIICF_P7_r3))
data[row.names(data) == "sp|P08047|SP1_HUMAN",]$LFS_05F_24_pre_r3 <- mean(c(data[row.names(data) == "sp|P08047|SP1_HUMAN",]$LFS_05F_24_pre_r1, data[row.names(data) == "sp|P08047|SP1_HUMAN",]$LFS_05F_24_pre_r2))

stopifnot(!is.na(data[row.names(data) == "sp|Q9C0F1|CEP44_HUMAN",]$WI38_r1))
stopifnot(!is.na(data[row.names(data) == "sp|P40305|IFI27_HUMAN",]$IIICF_P7_r2))
stopifnot(!is.na(data[row.names(data) == "sp|P08047|SP1_HUMAN",]$LFS_05F_24_pre_r3))
stopifnot(sum(is.na(data)) == 0)

data <- log2(data) # Absolute values are quite large ...

tests <- getTests()
pData <- getPData()

# Qunatile normalise log2 transformed data
normalised_exprs_log <- as.matrix(normalizeQuantiles(data, ties=TRUE))

# Create an object that LIMMA needs as input
protein_level_data_expression <- ExpressionSet(assayData=normalised_exprs_log, phenoData=pData, annotation="human")

# Box plots before quantile normalistaion
boxplot(data, main="Before Quantile Normalisation")

# Box plots after quantile normalistaion
boxplot(exprs(protein_level_data_expression), main="After Quantile Normalisation")

diff.test <- NULL
for (i in 1:nrow(tests))
{
    print(i)
    r <- runTest(data, tests[i,])
    if (is.null(r)) print(tests[i,]) else diff.test <- rbind(diff.test, r)
}

diff.test <- diff.test[with(diff.test, order(qval)),]
write.table(data.norm, file="SWATH/diffTest.tsv", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
