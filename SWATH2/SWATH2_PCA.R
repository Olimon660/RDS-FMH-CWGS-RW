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

   # if (showShape) {
        g <- autoplot(prcomp(data[,1:nc]), data=data, loadings.colour="blue", loadings.label=FALSE, loadings.label.size=0.1, size=1, label=FALSE, shape=TRUE)
 #   } else {
#        g <- autoplot(prcomp(data[,1:nc], na.action = na.omit), data=data, loadings.colour="blue", loadings.label=FALSE, loadings.label.size=0.1, size=1)
  #  }
    
    g <- g + geom_point(size=2.5, aes(col=Cell, shape=factor(Mortality)))
    g <- g + theme_bw()
    g <- g + ggtitle(title) +  theme(legend.title=element_blank())
    
    if (!showLeg) { g <- g + theme(legend.position="none") }
    print(g)
}

data <- read.table("SWATH2/SWATH2_data.tsv", header=TRUE, sep='\t', check.names=F)
data <- data[rowSums(is.na(data[,7:ncol(data)])) != (ncol(data) - 7),] # Remove all missing values
data <- log2(data[,7:ncol(data)])
info <- read.table("SWATH2/SWATH2_track.tsv", header=TRUE, sep='\t')
data[is.na(data)] <- 0

plotPCA(data, info, "Log2 PCA before normalization (colored by cells) (Peptide level)") # Before normalization








#normalize <- function(data)
#{
#    x <- normalize.quantiles(as.matrix(data))
#    colnames(x) <- colnames(data)
#    as.data.frame(x)
#}

#data.norm <- normalize(data)
#plotPCA(data.norm, info, "Log2 PCA after normalization (colored by cells)", "Cell")

#write.table(data, file="SWATH/beforeData.tsv", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
#write.table(data.norm, file="SWATH/afterData.tsv", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
