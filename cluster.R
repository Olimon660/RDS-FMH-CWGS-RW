library(VariantAnnotation)

path <- '/Users/twong/Sources/RDS-FMH-CWGS-RW/6'

parse <- function(x)
{
    vcf <- readVcf(x, "hg19")    
    CHR <- as.character(seqnames(vcf))
    POS <- start(ranges(vcf))
    REF <- as.character(ref(vcf))
    ALT <- as.character(unlist(alt(vcf)))
    key <- paste(CHR,POS,REF,ALT,sep='_')
    key
}

setwd(path)
files <- list.files(path, pattern = "\\vcf$")
data <- lapply(files, parse)

d <- proxy::dist(data, method = "Jaccard")
h = hclust(d)
plot(h)
