library(plyr)

getTests <- function()
{
    tests <- read.table("SWATH2/contrasts.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
    tests[tests$Mortal == "IIICF",]$Mortal <- "IIICF_P7" # IIICF_P7 is being used as the reference
    tests[tests$Mortal == "JFCF6_P12",]$Mortal <- "JFCF_6"
    tests    
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

trk <- merge(data.frame(IntenSample=colnames(inten)), trk, by.x="IntenSample", by.y="FileID")
write.table(trk$Sample, "/tmp/A.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
system("python3 SWATH2/SWATH2_COMPUTER.py > /tmp/B.txt")
trk$Sample <- read.table("/tmp/B.txt")$V1

trk_ <- subset(trk, grepl("JFCF", trk$Sample) | grepl("IIICF", trk$Sample) | grepl("GM847", trk$Sample) | grepl("IVG", trk$Sample) |
                    grepl("LFS", trk$Sample)  | grepl("MeT", trk$Sample)   | grepl("VA13", trk$Sample)  | grepl("GM02063", trk$Sample) |
                    grepl("WI38", trk$Sample))
inten_ <- inten[, colnames(inten) %in% trk_$ID]
inten_ <- cbind(info, inten_)

write.table(trk_, file="SWATH2/SWATH2_track.tsv", row.names=FALSE, quote=FALSE, sep='\t')
write.table(inten_, file="SWATH2/SWATH2_data.tsv", row.names=FALSE, quote=FALSE, sep='\t')
write.table(getTests(), file="SWATH2/SWATH2_tests.tsv", row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE)
