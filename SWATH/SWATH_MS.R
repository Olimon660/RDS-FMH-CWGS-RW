library(MSstats)

data <- read.csv("SWATH/Observed RT.csv")
saveRDS(data, file="SWATH/ObservedRT.rds")
data <- readRDS("SWATH/ObservedRT.rds")

#
# We need "FragmentIon", "sIsotopeLabelType",
#         "Condition", "BioReplicate", "Run", "Intensity".
#

track <- track()
info  <- data[,c(1:7)]
inten <- data[,c(8:ncol(data))]

info$ProductCharge <- NA
names(info)[names(info) == "Protein"] <- "ProteinName"
names(info)[names(info) == "Peptide"] <- "PeptideSequence"
names(info)[names(info) == "Precursor.Charge"] <- "PrecursorCharge"

colnames(inten) <- computerFriendly(colnames(inten))

#colnames(inten) %in% track$SWATH.File.Name

