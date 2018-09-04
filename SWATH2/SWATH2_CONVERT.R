library(MSstats)

data <- read.csv("SWATH2/Cell_Study_Reanalysis_All_415_samples_August_2018_PGH_peptide_FDR.csv", header=1)



#
# We need "FragmentIon", "IsotopeLabelType",
#         "Condition", "BioReplicate", "Run"
#
#  IsotopeLabelType not used
#

trk <- track()
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
stopifnot(sum(!(colnames(inten) %in% trk$FileID)) == 0)
