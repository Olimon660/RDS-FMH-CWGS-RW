library(plyr)

getTests <- function()
{
    test <- read.table("SWATH2/contrasts.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
    test[test$Mortal == "IIICF",]$Mortal <- "IIICF_P7" # IIICF_P7 is being used as the reference
    test[test$Mortal == "JFCF6_P12",]$Mortal <- "JFCF_6"
    colnames(test) <- c("Mortal", "Immortal")    
    test
}

track <- function()
{
    trk <- read.table("SWATH2/Cell Study_Tracking_SM.csv", header=FALSE, sep=',', stringsAsFactors=FALSE)
    trk <- trk[,c(1:15)]
    colnames(trk) <- c("ID", "Sample", "ProcessDate", "AcqType", "ProcessInst", "FileName", "Location", "Notes", "RunDate", "IDAFile", "IDALocation",
                       "MSUsed", "Operator", "MSMethod", "AcqType")
    trk <- trk[, c("ID", "Sample", "ProcessDate", "AcqType", "ProcessInst", "FileName", "RunDate", "MSUsed", "Operator", "MSMethod")]
    trk <- trk[with(trk, order(Sample, FileName)),]
    trk <- trk[trk$ID != "Sample No ",]
    trk <- trk[trk$AcqType != "IDA",]
    trk <- trk[trk$FileName != "",]
    trk$Sample <- trimws(trk$Sample)
    trk
}

data <- read.csv("SWATH2/Cell_Study_Reanalysis_All_415_samples_August_2018_PGH_peptide_FDR.csv", header=1)

trk <- track()
write.table(trk$FileName, "/tmp/A.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
system("python3 SWATH2/SWATH2_SAMPLE.py /tmp/A.txt /tmp/B.txt")
trk$FileID <- read.table("/tmp/B.txt")$V1

info  <- data[,c(1:6)]
inten <- data[,c(7:ncol(data))]
write.table(colnames(inten), "/tmp/A.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
system("python3 SWATH2/SWATH2_SAMPLE.py /tmp/A.txt /tmp/B.txt")
colnames(inten) <- read.table("/tmp/B.txt")$V1
stopifnot(sum(!(colnames(inten) %in% trk$FileID)) == 0)

tmp <- merge(data.frame(IntenSample=colnames(inten)), trk, by.x="IntenSample", by.y="FileID")
write.table(tmp$Sample, "/tmp/A.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
system("python3 SWATH2/SWATH2_COMPUTER.py > /tmp/B.txt")
tmp$FriendSample <- read.table("/tmp/B.txt", stringsAsFactors=F)$V1
tmp <- tmp[tmp$FriendSample != "????",]

trk_ <- subset(tmp, grepl("JFCF", tmp$FriendSample) | grepl("IIICF", tmp$FriendSample) | grepl("GM847", tmp$FriendSample) | grepl("IVG", tmp$FriendSample) |
                    grepl("LFS", tmp$FriendSample)  | grepl("MeT", tmp$FriendSample)   | grepl("VA13", tmp$FriendSample)  | grepl("GM02063", tmp$FriendSample) |
                    grepl("WI38", tmp$FriendSample))
trk_$Cell <- sapply(strsplit(trk_$FriendSample, "_"), `[`, 1)

inten_ <- inten[, colnames(inten) %in% trk_$ID]
inten_ <- cbind(info, inten_)

mortals <- c("JFCF_6", "GM02063", "IIICF_E6E7_C4_pre", "IVG_BF_LXSN_pre", "LFS_05F_24_pre", "MeT_4A_pre", "WI38", "IIICF_P7", "IIICF_P9")
trk_$Mortality <- "Immortal"
trk_[trk_$FriendSample %in% mortals,]$Mortality <- "Mortal"

write.table(trk_, file="SWATH2/SWATH2_track.tsv", row.names=FALSE, quote=FALSE, sep='\t')
write.table(inten_, file="SWATH2/SWATH2_data.tsv", row.names=FALSE, quote=FALSE, sep='\t')
write.table(getTests(), file="SWATH2/SWATH2_tests.tsv", row.names=FALSE, quote=FALSE, sep='\t')
