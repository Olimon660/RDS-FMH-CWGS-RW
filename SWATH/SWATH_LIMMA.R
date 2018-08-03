library(limma)
library(qvalue)
library(Biobase)
source("SWATH/SWATH_LIMMA_TOOLS.R")

data <- read.table("SWATH/data2.tsv", row.names=1, header=TRUE, sep='\t')
data <- data[!is.na(rowSums(data)),] # Should we do imputation?
data <- log2(data) # Absolute values are quite large ...

noRep <- function(x) { gsub('_r1', '', gsub('_r2', '', gsub('_r3', '', x))) }

# Qunatile normalise log2 transformed data
normalised_exprs_log <- as.matrix(normalizeQuantiles(data, ties=TRUE))

pData <- data.frame(row.names=colnames(data), Cond=noRep(colnames(data)))

# Create a AnnotatedDataFrame object
phenoData <- new("AnnotatedDataFrame", data=pData)

# Create an object that LIMMA needs as input
protein_level_data_expression <- ExpressionSet(assayData=normalised_exprs_log, phenoData=phenoData, annotation="human")

# Box plots before quantile normalistaion
boxplot(data, main="Before Quantile Normalisation")

# Box plots after quantile normalistaion
boxplot(exprs(protein_level_data_expression), main="After Quantile Normalisation")

# print the dimensions of the expression object form our main working data set 
dim(exprs(protein_level_data_expression))

tests <- read.table("SWATH/contrasts.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
tests[tests$Mortal == "IIICF",]$Mortal <- "IIICF_P7" # IIICF_P7 is being used as the reference

runTest <- function(data, test)
{
    case <- data[, noRep(colnames(data)) == test$Mortal[[1]]]   # Mortal
    cont <- data[, noRep(colnames(data)) == test$Immortal[[1]]] # Immortal

    if (ncol(case) != 3) { return(NULL) }
    if (ncol(cont) != 3) { return(NULL) }    
    
    stopifnot(length(unique(noRep(colnames(case)))) == 1)
    stopifnot(length(unique(noRep(colnames(cont)))) == 1)    
    
    # Data that used for testing
    newData <- cbind(case, cont)
    
    # http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
    design <- model.matrix(~factor(c(2,2,2,1,1,1)))
    colnames(design) <- c("Intercept", "Diff")
    
    eb.fit(newData, design)
}

res <- NULL
for (i in 1:nrow(tests))
{
    r <- runTest(data, tests[i,])
    if (is.null(r)) print(tests[i,]) else res <- rbind(res, r)
}
