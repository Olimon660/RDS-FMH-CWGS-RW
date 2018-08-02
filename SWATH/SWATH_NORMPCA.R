library(ggfortify)
library(data.table)

normalize <- function(data)
{
    x <- normalize.quantiles(as.matrix(data[,c(5:22)]))
    x <- cbind(data[,c(1:4)], data.frame(x))
    colnames(x) <- colnames(data)
    x
}

samples <- function(x)
{
    x <- gsub('_r1', '',  x)
    x <- gsub('_r2', '',  x)
    x <- gsub('_r3', '',  x)    
    x <- gsub('_r4', '',  x)    
    x    
}

plotPCA <- function(data, info, title, mode)
{ 
    rn <- colnames(data)
    data <- transpose(data)
    nc <- ncol(data)
    
    data$Samp <- rn
    data$Cond <- samples(rn) # Sample names without replicates
    stopifnot(all(sort(data$Samp) == sort(info$Sample)))
    
    # Merge with information such as instruments
    info <- merge(data, info, by.x="Samp", by.y="Sample")

    data$Inst <- info$Instrument # Useful for plotting
    
    g <- autoplot(prcomp(data[,1:nc]), data=data, loadings.colour="blue", loadings.label=FALSE, loadings.label.size=0.1, size=1.0)
    if (mode == "Samples")     { g <- g + geom_point(aes(col=Samp)) }
    if (mode == "Instruments") { g <- g + geom_point(aes(col=Inst)) }
    g <- g + theme_bw()
    g <- g + ggtitle(title) + guides(fill=FALSE) + theme(legend.title=element_blank())
    print(g)
}

data <- read.table("/Users/twong/Sources/RDS-FMH-CWGS-RW/SWATH/data2.tsv", row.names=1, header=TRUE, sep='\t')
info <- read.table("/Users/twong/Sources/RDS-FMH-CWGS-RW/SWATH/newSWATHDetails.tsv", header=TRUE, sep='\t')

data <- data[!is.na(rowSums(data)),] # Should we do imputation?
data <- log2(data) # The absolute values are quite large ...

plotPCA(data, info, "PCA before normalization (colored by samples)", "Samples")
plotPCA(data, info, "PCA before normalization (colored by instruments)", "Instruments")
