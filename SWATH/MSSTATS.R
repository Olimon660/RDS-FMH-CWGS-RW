library(MSstats)

data <- read.csv("SWATH/Observed RT.csv")
saveRDS(data, file="SWATH/ObservedRT.rds")
data <- readRDS("SWATH/ObservedRT.rds")

#
# We need "FragmentIon", "sIsotopeLabelType",
#         "Condition", "BioReplicate", "Run", "Intensity".
#

data$ProductCharge <- NA
names(data)[names(data) == "Protein"] <- "ProteinName"
names(data)[names(data) == "Peptide"] <- "PeptideSequence"
names(data)[names(data) == "Precursor.Charge"] <- "PrecursorCharge"

dataProcess(data)
