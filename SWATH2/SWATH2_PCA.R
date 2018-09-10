library(ggfortify)
library(data.table)
library(preprocessCore)

plotPCA <- function(data, info, title, showShape=FALSE, showLeg=TRUE)
{ 
    cn <- colnames(data)
    data <- transpose(data)
    nc <- ncol(data)
    
    data$ID <- cn
    data_ <- merge(data, info, by.x="ID", by.y="ID")
    stopifnot(nrow(data) == nrow(data_))
    data <- data_
    data$ID <- NULL

    g <- autoplot(prcomp(data[,1:nc]), data=data, loadings.colour="blue", loadings.label=FALSE, loadings.label.size=0.1, size=1, label=FALSE, shape=TRUE)
    g <- g + geom_point(size=2.5, aes(col=Cell, shape=factor(Mortality)))
    g <- g + theme_bw()
    g <- g + ggtitle(title) +  theme(legend.title=element_blank())
    
    if (!showLeg) { g <- g + theme(legend.position="none") }
    print(g)
}

toProts <- function(pd, info)
{
    aggregate(as.matrix(pd[,7:ncol(pd)]) ~ Protein, data=pd, sum, na.rm=TRUE)
}

data <- read.table("SWATH2/SWATH2_data.tsv", header=TRUE, sep='\t', check.names=F)
data <- data[rowSums(is.na(data[,7:ncol(data)])) != (ncol(data) - 7),] # Remove all missing values
data[,7:ncol(data)] <- log2(data[,7:ncol(data)])
info <- read.table("SWATH2/SWATH2_track.tsv", header=TRUE, sep='\t')

pd <- data         # data at the peptide level
pd[is.na(pd)] <- 0 # Also useful for the next step (no missing value)
pt <- toProts(pd, info)

plotPCA(pd[,7:ncol(pd)], info, "Log2 PCA before normalization (colored by cells) (Peptide level)") # Before normalization
plotPCA(pt[,2:ncol(pt)], info, "Log2 PCA before normalization (colored by cells) (Protein level)") # Before normalization

normalize <- function(data)
{
    x <- normalize.quantiles(as.matrix(data))
    colnames(x) <- colnames(data)
    as.data.frame(x)
}

pd.norm <- pd
pt.norm <- pt
pd.norm[,7:ncol(pd.norm)] <- normalize(pd.norm[,7:ncol(pd.norm)])
pt.norm[,2:ncol(pt.norm)] <- normalize(pt.norm[,2:ncol(pt.norm)])

plotPCA(pd.norm[,7:ncol(pd.norm)], info, "Log2 PCA after normalization (colored by cells) (Peptide level)") # After normalization
plotPCA(pt.norm[,2:ncol(pt.norm)], info, "Log2 PCA after normalization (colored by cells) (Protein level)") # After normalization
