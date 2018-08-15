library(MSstats)

#data <- read.csv("SWATH/Observed RT.csv")
#saveRDS(data, file="SWATH/ObservedRT.rds")
#data <- readRDS("SWATH/ObservedRT.rds")

#
# We need "FragmentIon", "sIsotopeLabelType",
#         "Condition", "BioReplicate", "Run", "Intensity".
#

trk <- track()
trk <- trk[trk$AcqType != "IDA",]
write.table(trk$FileName, "/tmp/A.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
system("python3 SWATH/SWATH_SAMPLE.py /tmp/A.txt /tmp/B.txt")
trk$FileID <- read.table("/tmp/B.txt")$V1

info  <- data[,c(1:7)]
inten <- data[,c(8:ncol(data))]

info$ProductCharge <- NA
names(info)[names(info) == "Protein"] <- "ProteinName"
names(info)[names(info) == "Peptide"] <- "PeptideSequence"
names(info)[names(info) == "Precursor.Charge"] <- "PrecursorCharge"

write.table(colnames(inten), "/tmp/A.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
system("python3 SWATH/SWATH_SAMPLE.py /tmp/A.txt /tmp/B.txt")
colnames(inten) <- read.table("/tmp/B.txt")$V1




    print(colnames(inten)[!(colnames(inten) %in% trk$FileID)])
#print(trk$FileName)


