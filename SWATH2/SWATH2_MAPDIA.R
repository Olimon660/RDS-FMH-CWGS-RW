library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)

info <- read.table("SWATH2/SWATH2_track.tsv", header=1, sep='\t', stringsAsFactors=F)
data <- read.table("SWATH2/SWATH2_data.tsv", header=1, sep='\t', check.names=F)
test <- read.table("SWATH2/SWATH2_tests.tsv", header=1, sep='\t', stringsAsFactors=F)
data <- dplyr::select(data, -c("Unique_code", "Precursor.MZ", "Precursor.Charge", "RT"))
data <- plyr::rename(data, c('Protein'='ProteinName'))
data <- plyr::rename(data, c('Peptide'='PeptideSequence'))
data <- add_column(data, FragmentIon='', .after = 2)

plotVol <- function(path, mortal, immortal)
{
    saveTo <- paste(mortal, "-", immortal, ".png", sep="")
    
    out <- read.table(paste(path, "/analysis_output.txt", sep=""), sep="\t", header=1)
    out$color <- "black"
    if (nrow(out[out$FDR <= 0.05,]) > 0)  { out[out$FDR <= 0.05,]$color <- "red" }
    if (nrow(out[abs(out$log2FC) >= 3,]) > 0) { out[abs(out$log2FC) >= 3,]$color <- "orange" }
    if (nrow(out[out$FDR <= 0.05 & abs(out$log2FC) >= 3,]) > 0) { out[out$FDR <= 0.05 & abs(out$log2FC) >= 3,]$color <- "green" }
    out$color <- as.factor(out$color)    
    
    out$Name <- as.character(out$Protein)
    out$Name <- sapply(strsplit(out$Name, "|", fixed=TRUE), `[`, 3)

    png(paste("~/Desktop/", saveTo, sep=""), width=480*2)    
    g <- ggplot(out, aes(x=log2FC, y=-log2(FDR)))
    g <- g + geom_point(aes(col=color))
    g <- g + scale_colour_manual(values = c(levels(out$color)))
    g <- g + geom_vline(xintercept=3, color="black", size=0.2)
    g <- g + geom_vline(xintercept=-3, color="black", size=0.2)
    g <- g + theme_bw()
    g <- g + xlab("Log fold changes")
    g <- g + ylab("-log2 FDR")
    g <- g + geom_text(data=out[out$color == "green",], aes(y=-log2(FDR), x=log2FC, label=Name), size=3.0)
    g <- g + ggtitle(paste("Mortal: ", mortal, "    Immortal: ", immortal, sep=""))
    g <- g + theme(legend.position="none")
    print(g)
    dev.off()
}

r <- NULL

for (row in 1:nrow(test))
{
    t <- test[row,]
    m <- info[info$FriendSample == t$Mortal,]   # Mortal samples
    i <- info[info$FriendSample == t$Immortal,] # Immortal samples

    if (t$Immortal == "JFCF_6_P_pLKO_5") {
        next
    }

    stopifnot(nrow(m) > 0)
    stopifnot(nrow(i) > 0)

    d <- data[, colnames(data) %in% m$ID | colnames(data) %in% i$ID | colnames(data) %in% c("ProteinName", "PeptideSequence", "FragmentIon")]
    colnames(d) <- c(c("ProteinName", "PeptideSequence", "FragmentIon"), c(rep(t$Mortal, nrow(m))), rep(t$Immortal, nrow(i)))

    n1 <- nrow(m) + nrow(i) # Total number of samples
    
    # Remove peptides with all missing data
    d <- d[rowSums(is.na(d[,c(4:ncol(d))])) != n1,]
    
    write.table(d, "/tmp/data.txt", quote=FALSE, row.names=FALSE, sep='\t')
    system("cp SWATH2/SWATH2_MAPDIA.txt /tmp/")
    system("rm analysis_output.txt; rm analysis_output_wide_format.txt")
    system("rm fragment_selection.txt; rm log2_data.txt")
    system("rm param.txt; rm protein_level.txt")
    system("mapDIA /tmp/SWATH2_MAPDIA.txt")

    t2 <- read.table(paste(getwd(), "/analysis_output.txt", sep=""), sep="\t", header=1)
    t2$Mortal <- t$Mortal
    t2$Immortal <- t$Immortal
    
    r <- rbind(r, data.frame(Mortal=t2$Mortal,
                             Immortal=t2$Immortal,
                             Protein=t2$Protein,
                             nPeptide=t2$nPeptide,
                             log2FC=t2$log2FC,
                             log2FC_SE=t2$log2FC_SE,
                             FDR=t2$FDR,
                             log_oddsDE=t2$log_oddsDE))

   plotVol(getwd(), t$Mortal, t$Immortal)
}



