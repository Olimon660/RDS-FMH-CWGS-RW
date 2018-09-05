library(plyr)
library(dplyr)
library(tibble)

#info <- read.table("SWATH2/SWATH2_track.tsv", header=1, sep='\t', stringsAsFactors=F)
#data <- read.table("SWATH2/SWATH2_data.tsv", header=1, sep='\t', check.names=F)
#test <- read.table("SWATH2/SWATH2_tests.tsv", header=1, sep='\t', stringsAsFactors=F)
#data <- dplyr::select(data, -c("Unique_code", "Precursor.MZ", "Precursor.Charge", "RT"))
#data <- plyr::rename(data, c('Protein'='ProteinName'))
#data <- plyr::rename(data, c('Peptide'='PeptideSequence'))
#data <- add_column(data, FragmentIon='', .after = 2)

for (row in 1:nrow(test))
{
    t <- test[row,]
    m <- info[info$Sample == t$Mortal,]   # Mortal samples
    i <- info[info$Sample == t$Immortal,] # Immortal samples

    stopifnot(nrow(m) > 0)
    stopifnot(nrow(i) > 0)

    d <- data[, colnames(data) %in% m$ID | colnames(data) %in% i$ID | colnames(data) %in% c("ProteinName", "PeptideSequence", "FragmentIon")]
    colnames(d) <- c(c("ProteinName", "PeptideSequence", "FragmentIon"), c(rep(t$Mortal, nrow(m))), rep(t$Immortal, nrow(i)))

    n1 <- nrow(m) + nrow(i) # Total number of samples
    
    # Remove peptides with all missing data
    d <- d[rowSums(is.na(d[,c(4:9)])) != n1,]
    
    # Write data file for mapDIA differential testing    
    write.table(d, "/tmp/data.txt", quote=FALSE, row.names=FALSE, sep='\t')

    break
}
