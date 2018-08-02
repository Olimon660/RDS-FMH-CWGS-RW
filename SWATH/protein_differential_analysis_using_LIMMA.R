# Authors: Erdahl Teber and Matloob Khushi
# Date: 13 June 2014
# Last modified: 13 June 2014
# Purpose: To perform protein level differential expression analysis for Peter Haines.


#reading in raw peak areas

file <- '/mnt/san/bioinformatics/mkhushi/PeterHains/Hains_RAW.txt'

tb <- read.table(file, header = TRUE, sep = "\t", quote = "\"",na.strings = "", as.is=TRUE)

# show firt 10 observations
head(tb)


library(sqldf)
# collapse data to protein level by taking the median value from all associated peptides 
# for each of the experiments (16 experiments).

protein_level <-sqldf("select Protein, median(SC2) as SC2, median(SC3) as SC3,
               median(SC4) as SC4, median(SC5) as SC5, median(SC6) as SC6,
               median(SD2) as SD2, median(SD3) as SD3,
               median(SD4) as SD4, median(SD5) as SD5, median(SD6) as SD6,
               
               median(SW1) as SW1, median(SW2) as SW2, median(SW3) as SW3,
               median(SW4) as SW4, median(SW5) as SW5, median(SW6) as SW6
               from tb group by Protein order by Protein")

# create a matrix object cntaining experimental peak areas 
exprs <- as.matrix(protein_level[ ,2:17], as.is=TRUE)

# log 2 transform 
exprs_log<-log2(exprs)
# qunatile normalise log2 transformed data
normalised_exprs_log <-normalizeQuantiles(exprs_log, ties=TRUE)

# assign protein ID to matrix. There each row is a unique protein ID
rownames(normalised_exprs_log)<- protein_level[,1]

#print column names - experimental columns of normalised protein level data
colnames(normalised_exprs_log)
head(normalised_exprs_log)

library(Biobase)

# read in phenotype data
pDataFile <- '/mnt/san/bioinformatics/mkhushi/PeterHains/pData.txt'
# second column is used for rownames
pData <- read.table(pDataFile, row.names=2, header = TRUE, sep = "\t", quote = "\"")

# create a AnnotatedDataFrame object
phenoData <- new("AnnotatedDataFrame", data=pData)
# print phenoData
phenoData

# create an object that LIMMA needs as input

protein_level_data_expression <- ExpressionSet(assayData=normalised_exprs_log, phenoData=phenoData, annotation="human")

# print box plots before quantile normalistaion

boxplot(exprs_log,main="Before Quantile Normalisation")

# print box plots after quantile normalistaion

boxplot(exprs(protein_level_data_expression),main="After Quantile Normalisation")

# print the dimensions of the expression object form our main working data set 

dim(exprs(protein_level_data_expression))

library(limma)

# create a design matrix for comparing between different experiments
# sampleTypes (control, addicted, withdrawal)

treatment <- factor(pData(protein_level_data_expression)[,"sampleTypes"])

design <- model.matrix(~0+treatment)
colnames(design)<- levels(treatment)
colnames(design)
# take into account any subtle inconsistenties between experiments, 
#by giving stronger weight to more highly correlated experiments within a treatment type (control, addicted, withdrawal)
aw <- arrayWeights(exprs(protein_level_data_expression), design)
# print the weights
aw

# fit a linear model using our comparison matrix and column (experiments) weights
fit <- lmFit(exprs(protein_level_data_expression), design, weights=aw)

names(fit)
dim(fit)
# set up the contrasts
contrasts <- makeContrasts(Addicted - Control, Withdrawal - Control, Addicted - Withdrawal, levels = design)
# use LIMMA empirical Bayes statistical technique to smooth out the variation 
contr.fit    <- eBayes(contrasts.fit(fit, contrasts))

# list the top 100 diffrentially expressed proteinIDs, with a adjusted p-value < 0.05

Addicted_Control_results        <-topTable(contr.fit,coef='Addicted - Control', sort="p", n=100, p.value = 0.05, lfc=0)
Withdrawal_Control_results      <-topTable(contr.fit,coef='Withdrawal - Control', sort="p", n=100, p.value = 0.05, lfc=0)
Addicted_Withdrawal_results     <-topTable(contr.fit,coef='Addicted - Withdrawal', sort="p", n=100, p.value = 0.05, lfc=0)

# do a MDS plot of all experiemnts using the top 500 proteins with the highest level of variation (i.e wide peak area distrubitions).
plotMDS(main='peak area similarities', exprs(protein_level_data_expression), 
        top=500, gene.selection="common",
        labels=rownames(pData(protein_level_data_expression)), 
        col=c('red','blue','green')[factor(pData(protein_level_data_expression)[,'sampleTypes'])])


# produce a vennDiagram of results
results <- decideTests(contr.fit, method="separate", adjust.method="BH", p.value=0.05,lfc=0)
vennDiagram(results, include=c("up","down"),counts.col=c("red","green"))  
