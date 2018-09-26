
#temp  <- read.csv("/Users/twong/Sources/RDS-FMH-CWGS-RW/SWATH2/Cell_Study_Reanalysis_All_415_samples_August_2018_PGH_peptide_FDR.csv")
#ATRX_ <- temp[grepl("DAXX", temp$Protein),]
#ATRX_ <- ATRX_[,7:ncol(ATRX_)]

data <- read.table("/Users/twong/Sources/RDS-FMH-CWGS-RW/SWATH2/SWATH2_data.tsv", sep="\t", header=1)
info <- read.table("/Users/twong/Sources/RDS-FMH-CWGS-RW/SWATH2/SWATH2_track.tsv", sep="\t", header=1)
ATRX <- data[grepl("ATRX", data$Protein),]
DAXX <- data[grepl("DAXX", data$Protein),]
TP53 <- data[grepl("TP53", data$Protein),]

ATRX_ <- data.frame(Protein=rep("ATRX", nrow(ATRX)), Peptide=ATRX$Peptide)
DAXX_ <- data.frame(Protein=rep("DAXX", nrow(DAXX)), Peptide=DAXX$Peptide)
TP53_ <- data.frame(Protein=rep("TP53", nrow(TP53)), Peptide=TP53$Peptide)

write.table(rbind(ATRX_, DAXX_, TP53_), "/Users/twong/Desktop/A.tsv", sep="\t", quote=FALSE, row.names=FALSE)

sum(!is.na((ATRX[,7:ncol(ATRX)])))
sum(!is.na((DAXX[,7:ncol(DAXX)])))
sum(!is.na((TP53[,7:ncol(TP53)])))
quantile(ATRX[,7:ncol(ATRX)], na.rm=TRUE)
quantile(DAXX[,7:ncol(DAXX)], na.rm=TRUE)
quantile(TP53[,7:ncol(TP53)], na.rm=TRUE)

#print(nrow(ATRX)) # 9
#print(nrow(DAXX)) # 5
#print(nrow(TP53)) # 32

tests <- read.csv("/Users/twong/Sources/RDS-FMH-CWGS-RW/AN1/AN1_CONTRASTS.csv", stringsAsFactors=FALSE)
r <- NULL

for (i in 1:nrow(tests))
{
    t  <- tests[i,]
    m1 <- t$Mortal
    m2 <- t$Immortal
    if (m1 == "IIICF")     { m1 <- "IIICF_P7" }
    if (m2 == "IIICF")     { m2 <- "IIICF_P7" }
    if (m1 == "JFCF6_P12") { m1 <- "JFCF_6"   }
    if (m2 == "JFCF6_P12") { m2 <- "JFCF_6"   }
    
    i1 <- info[info$FriendSample == m1,]
    i2 <- info[info$FriendSample == m2,]
    print(m1)
    
    stopifnot(nrow(i1) > 0)
    stopifnot(nrow(i2) > 0)
    
    x1 <- paste("X", i1$ID, sep="")
    x2 <- paste("X", i2$ID, sep="")
    d1 <- ATRX[,colnames(ATRX) %in% x1]
    d2 <- ATRX[,colnames(ATRX) %in% x2]
    d1[is.na(d1)] <- 0
    d2[is.na(d2)] <- 0
    c1 <- colSums(d1)
    c2 <- colSums(d2)
    a1 <- mean(c1)
    a2 <- mean(c2)
    r  <- rbind(r, data.frame(A1=a1, A2=a2, Mortal=m1, Immortal=m2))       
}