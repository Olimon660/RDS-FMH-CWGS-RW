library(plyr)
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

trk  <- merge(data.frame(IntenSample=colnames(inten)), trk, by.x="IntenSample", by.y="FileID")
trk_ <- subset(trk, grepl("JFCF", trk$Sample) | grepl("IIICF", trk$Sample) | grepl("GM847", trk$Sample) | grepl("IVG", trk$Sample) |
                    grepl("LFS", trk$Sample)  | grepl("MeT", trk$Sample)   | grepl("VA13", trk$Sample)  | grepl("GM02063", trk$Sample) |
                    grepl("WI38", trk$Sample))
inten_ <- inten[, colnames(inten) %in% trk_$ID]

write.table(trk_$Sample, "/tmp/A.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
system("python3 SWATH2/SWATH2_COMPUTER.py > /tmp/B.txt")
trk_$Computer <- read.table("/tmp/B.txt")$V1

write.table(trk_,   file="SWATH2/SWATH2_track.tsv", header=FALSE)
write.table(inten_, file="SWATH2/SWATH2_intensity.tsv", header=FALSE)
