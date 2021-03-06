library(ggfortify)
library(data.table)
library(preprocessCore)

noRep <- function(x) { gsub('_r1', '', gsub('_r2', '', gsub('_r3', '', x))) }

plotPCA <- function(data, info, title, mode, showShape=FALSE, showLeg=TRUE)
{ 
    rn <- colnames(data)
    data <- transpose(data)
    nc <- ncol(data)
    
    data$Samp <- rn
    data$Cond <- noRep(rn) # Sample names without replicates
    stopifnot(all(sort(data$Samp) == sort(info$Sample)))
    
    # Merge with information such as instruments
    info <- merge(data, info, by.x="Samp", by.y="Sample")

    data$Inst <- info$Instrument
    data$RunDate <- info$RunDate
    data$ProcessDate <- info$ProcessDate
    data$Cell <- info$Cell
    data$Type <- info$Type
    data$MSUsed <- info$MSUsed
    data$MSMethod <- info$MSMethod
    data$Operator <- info$Operator
    data$Mortality <- info$Mortality
    
    if (showShape) {
        g <- autoplot(prcomp(data[,1:nc]), data=data, loadings.colour="blue", loadings.label=TRUE, loadings.label.size=0.1, size=1, label=TRUE, shape=TRUE)
    } else {
        g <- autoplot(prcomp(data[,1:nc]), data=data, loadings.colour="blue", loadings.label=FALSE, loadings.label.size=0.1, size=1)
    }

    if (mode == "Cell")        { g <- g + geom_point(size=2.5, aes(col=Cell, shape=factor(Mortality))) }
    if (mode == "Samples")     { g <- g + geom_point(size=2.5, aes(size=0.3, col=Samp)) }
    if (mode == "Instruments") { g <- g + geom_point(size=2.5, aes(size=0.3, col=Inst)) }
    if (mode == "RunDate")     { g <- g + geom_point(size=2.5, aes(size=0.3, col=RunDate)) }
    if (mode == "ProcessDate") { g <- g + geom_point(size=2.5, aes(size=0.3, col=ProcessDate)) }
    if (mode == "Type")        { g <- g + geom_point(size=2.5, aes(size=0.3, col=Type)) }
    if (mode == "MSUsed")      { g <- g + geom_point(size=2.5, aes(size=0.3, col=MSUsed)) }
    if (mode == "MSMethod")    { g <- g + geom_point(size=2.5, aes(size=0.3, col=MSMethod)) }
    if (mode == "Operator")    { g <- g + geom_point(size=2.5, aes(size=0.3, col=Operator)) }
    if (mode == "Mortality")   { g <- g + geom_point(size=2.5, aes(size=0.3, col=Mortality)) }
    if (mode == "Condition")   { g <- g + geom_point(size=2.5, aes(size=0.3, col=Cond)) }
    
    g <- g + theme_bw()
    g <- g + scale_colour_manual(values = c("red", "orange", "brown", "green", "blue", "black", "violet"))
    g <- g + ggtitle(title) +  theme(legend.title=element_blank())

    if (!showLeg) { g <- g + theme(legend.position="none") }
    print(g)
}

data <- read.table("SWATH/data2.tsv", row.names=1, header=TRUE, sep='\t')
data <- data[!is.na(rowSums(data)),] # Should we do imputation?
data <- log2(data) # The absolute values are quite large ...
info <- read.table("SWATH/nInfo.tsv", header=TRUE, sep='\t')

plotPCA(data, info, "Log2 PCA before normalization (colored by cells)", "Cell")

normalize <- function(data)
{
    x <- normalize.quantiles(as.matrix(data))
    colnames(x) <- colnames(data)
    as.data.frame(x)
}

data.norm <- normalize(data)
plotPCA(data.norm, info, "Log2 PCA after normalization (colored by cells)", "Cell")

write.table(data, file="SWATH/beforeData.tsv", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(data.norm, file="SWATH/afterData.tsv", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
