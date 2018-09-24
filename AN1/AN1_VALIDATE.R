
info <- read.table("SWATH2/SWATH2_track.tsv", header=TRUE, sep='\t')

DAXX <- df1[grepl("DAXX", df1$Unique_code),]
ATRX$Protein <- factor(ATRX$Protein)
DAXX$Protein <- factor(DAXX$Protein)


df1  <- read.csv("SWATH2/SWATH2_data.tsv", header=1, sep="\t")
ATRX <- df1[grepl("ATRX", df1$Unique_code),]
a <- colSums(ATRX[,7:ncol(ATRX)], na.rm=TRUE)




DAXX[,7:ncol(DAXX)] <- sum(DAXX[,7:ncol(DAXX)], na.rm=TRUE)
