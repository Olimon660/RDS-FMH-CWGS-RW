library(limma)
library(qvalue)
library(Biobase)
source("SWATH/SWATH_LIMMA_TOOLS.R")

noRep <- function(x) { gsub('_r1', '', gsub('_r2', '', gsub('_r3', '', x))) }

data <- read.table("SWATH/data2.tsv", row.names=1, header=TRUE, sep='\t')
#data <- data[!is.na(rowSums(data)),] # Should we do imputation?
#data <- log2(data) # Absolute values are quite large ...

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

# print the dimensions of the expression object form our main working data set 
dim(exprs(protein_level_data_expression))


res <- NULL
for (i in 1:nrow(tests))
{
    r <- runTest(data, tests[i,])
    if (is.null(r)) print(tests[i,]) else res <- rbind(res, r)
}
