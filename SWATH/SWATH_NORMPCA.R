library(ggfortify)
library(data.table)
library(preprocessCore)

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

    data$Inst <- info$Instrument
    data$RunDate <- info$RunDate
    data$ProcessDate <- info$ProcessDate
    data$Type <- info$Type
    data$MSUsed <- info$MSUsed
    data$MSMethod <- info$MSMethod
    data$Operator <- info$Operator
    
    g <- autoplot(prcomp(data[,1:nc]), data=data, loadings.colour="blue", loadings.label=FALSE, loadings.label.size=0.1, size=1.0)
    
    if (mode == "Samples")     { g <- g + geom_point(aes(col=Samp)) }
    if (mode == "Instruments") { g <- g + geom_point(aes(col=Inst)) }
    if (mode == "RunDate")     { g <- g + geom_point(aes(col=RunDate)) }
    if (mode == "ProcessDate") { g <- g + geom_point(aes(col=ProcessDate)) }
    if (mode == "Type")        { g <- g + geom_point(aes(col=Type)) }
    if (mode == "MSUsed")      { g <- g + geom_point(aes(col=MSUsed)) }
    if (mode == "MSMethod")    { g <- g + geom_point(aes(col=MSMethod)) }
    if (mode == "Operator")    { g <- g + geom_point(aes(col=Operator)) }
    
    g <- g + theme_bw() + theme(legend.position="none")
    g <- g + ggtitle(title) + theme(legend.title=element_blank())
    print(g)
}

data <- read.table("SWATH/data2.tsv", row.names=1, header=TRUE, sep='\t')
info <- read.table("SWATH/nInfo.tsv", header=TRUE, sep='\t')
data <- data[!is.na(rowSums(data)),] # Should we do imputation?
data <- log2(data) # The absolute values are quite large ...

png("~/Desktop/BSamples.png");     plotPCA(data, info, "Log2 PCA before normalization (colored by samples)", "Samples");             dev.off()
png("~/Desktop/BInstruments.png"); plotPCA(data, info, "Log2 PCA before normalization (colored by instruments)", "Instruments");     dev.off()
png("~/Desktop/BProcessDate.png"); plotPCA(data, info, "Log2 PCA before normalization (colored by processing date)", "ProcessDate"); dev.off()
png("~/Desktop/BOperator.png");    plotPCA(data, info, "Log2 PCA before normalization (colored by operator)", "Operator");           dev.off()
png("~/Desktop/BRunDate.png");     plotPCA(data, info, "Log2 PCA before normalization (colored by running date)", "RunDate");        dev.off()
png("~/Desktop/BType.png");        plotPCA(data, info, "Log2 PCA before normalization (colored by type)", "Type");                   dev.off()
png("~/Desktop/BMSUsed.png");      plotPCA(data, info, "Log2 PCA before normalization (colored by MSUsed)", "MSUsed");               dev.off()
png("~/Desktop/BMSMethod.png");    plotPCA(data, info, "Log2 PCA before normalization (colored by MSMethod)", "MSMethod");           dev.off()

normalize <- function(data)
{
    x <- normalize.quantiles(as.matrix(data))
    colnames(x) <- colnames(data)
    as.data.frame(x)
}

data.norm <- normalize(data)

png("~/Desktop/ASamples.png");     plotPCA(data.norm, info, "Log2 PCA after normalization (colored by samples)", "Samples");             dev.off()
png("~/Desktop/AInstruments.png"); plotPCA(data.norm, info, "Log2 PCA after normalization (colored by instruments)", "Instruments");     dev.off()
png("~/Desktop/AProcessDate.png"); plotPCA(data.norm, info, "Log2 PCA after normalization (colored by processing date)", "ProcessDate"); dev.off()
png("~/Desktop/AOperator.png");    plotPCA(data.norm, info, "Log2 PCA after normalization (colored by operator)", "Operator");           dev.off()
png("~/Desktop/ARunDate.png");     plotPCA(data.norm, info, "Log2 PCA after normalization (colored by running date)", "RunDate");        dev.off()
png("~/Desktop/AType.png");        plotPCA(data.norm, info, "Log2 PCA after normalization (colored by type)", "Type");                   dev.off()
png("~/Desktop/AMSUsed.png");      plotPCA(data.norm, info, "Log2 PCA after normalization (colored by MSUsed)", "MSUsed");               dev.off()
png("~/Desktop/AMSMethod.png");    plotPCA(data.norm, info, "Log2 PCA after normalization (colored by MSMethod)", "MSMethod");           dev.off()
