
data <- read.table("/Users/twong/Sources/RDS-FMH-CWGS-RW/AN1/AN1_SUMMARY.tsv", sep="\t", header=1, stringsAsFactors=FALSE)

data$MStatus <- ifelse(data$GM_SNP > 0 | data$GM_Ind | data$VM_MUT > 0, "+", "-") # Mortal
data$IStatus <- ifelse(data$GI_SNP > 0 | data$GI_Ind | data$VI_MUT > 0, "+", "-") # Immortal

data$MStatus <- paste(data$Mortal,   data$MStatus, sep="")
data$IStatus <- paste(data$Immortal, data$IStatus, sep="")

write.table(data.frame(Mortal=data$MStatus, Immortal=data$IStatus), "/Users/twong/Sources/RDS-FMH-CWGS-RW/AN1/AN1_PATHWAY.tsv", sep="\t", quote=FALSE, row.names=FALSE)
