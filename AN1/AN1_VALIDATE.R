
data <- read.table("AN1/AN1_RESULTS.tsv", sep="\t", header=TRUE)

# No mutation for any gene/promoter in IIICF?
sum(data[data$Mortal == "IIICF",]$WGS_Mortal_SNP + data[data$Mortal == "IIICF",]$WGS_Mortal_Ind)
