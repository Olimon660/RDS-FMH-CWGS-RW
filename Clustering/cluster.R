library(circlize)
library(ComplexHeatmap)
library(VariantAnnotation)

#path <- '/Users/twong/Desktop/6'
path <- '/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/6'

setwd(path)
files <- list.files(path, pattern = "^DECOMPOSED_SNP_GATK_.*vcf$") # Only SNPs
data <- data.frame()

for (file in files)
{
    print(file)
    vcf <- readVcf(file, 'hg19')
    CHR <- as.character(seqnames(vcf))
    POS <- start(ranges(vcf))
    REF <- as.character(ref(vcf))
    ALT <- as.character(unlist(alt(vcf))) # Assume no multi-allele (decomposed)
    key <- paste(CHR, POS, REF, ALT, sep='_')
    data  <- rbind(data, data.frame(Key=key, Sample=file))
}

saveRDS(data, 'data.rds')
data <- readRDS('data.rds')

uniq <- unique(sort(data$Key))
print(paste(c('There are'), length(uniq), 'unique variants'))

m <- matrix(data=0, nrow=length(files), ncol=length(uniq)) # Samples vs SNPs
rownames(m) <- files
colnames(m) <- uniq

# Set 1 if the SNP appears in the sample (otherwise the matrix is 0)
ix <- unique(cbind(as.character(data$Sample), as.character(data$Key)))
m[ix] <- 1

saveRDS(m, 'm.rds')
m <- readRDS('m.rds')

jaccard <- function(M, samp1, samp2) {
    sums = colSums(M[c(samp1, samp2),])
    similarity = length(sums[sums==2])
    total = length(sums[sums==1]) + similarity
    similarity/total
}

jaccard_index<-matrix(data=0, nrow=length(files), ncol=length(files))
rownames(jaccard_index) <- files
colnames(jaccard_index) <- files

for (i in rownames(jaccard_index))
{
    for (j in rownames(jaccard_index))
    {
        n <- n + 1
        jaccard_index[i,j] <- jaccard(m,i,j)
    }
}

saveRDS(jaccard_index, 'jaccard_index.rds')
jaccard_index <- readRDS('jaccard_index.rds')

png('hclust.png')
plot(hclust(dist(jaccard_index)), xlab='')
dev.off()

fixNames <- function(x)
{
    x <- gsub('.vcf', '', gsub('GATK_', '', x))
    x <- gsub('_H06L4ALXX_3', '', x)
    x <- gsub('_H06L4ALXX_6', '', x)
    x <- gsub('_H06L4ALXX_4', '', x)
    x
}

col_text <- function(x) {
    if (x == 0) { y<-'' }
    else        { y<-x }
    return(y)
}

drawHeat <- function(m)
{
    Heatmap(m, name = "Jaccard index", cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 2, type = "none"), 
            col=colorRamp2(c(min(m), max(m)), c('white','red')), show_heatmap_legend = TRUE, 
            column_title = "SNV Indel similarity across immortal samples",
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width*0.90, height = height*0.90,gp = gpar(fill = fill))
                grid.text(col_text(round(m[i, j],2)), vjust=1, hjust=0.50, x, y, gp = gpar(fontsize = 8, col='black'))
            })
}

fixMatrix <- function(m, x)
{
    colnames(m) <- fixNames(colnames(m))
    rownames(m) <- fixNames(rownames(m))
    m <- m[!rownames(m) %in% x,]
    m <- m[,!colnames(m) %in% x]    
    m
}

drawHeat(fixMatrix(as.matrix(jaccard_index), c('ZK-58', 'VA13', 'KPD')))
