library(proxy)
library(VariantAnnotation)

#path <- '/Users/twong/Sources/RDS-FMH-CWGS-RW/6'
path <- '/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/6'

setwd(path)
files <- list.files(path, pattern = "\\vcf$")
data <- data.frame()

for (file in files)
{
    print(file)
    vcf <- readVcf(file, "hg19")
    CHR <- as.character(seqnames(vcf))
    POS <- start(ranges(vcf))
    REF <- as.character(ref(vcf))
    ALT <- as.character(unlist(alt(vcf)))
    key <- paste(CHR, POS, REF, ALT, sep='_')
    data  <- rbind(data, data.frame(Key=key, Sample=file))
}

uniq <- unique(sort(data$Key))
print(paste(c('There are'), length(uniq), 'unique variants'))

m <- matrix(data=0, nrow=length(files), ncol=length(uniq))
rownames(m) <- files
colnames(m) <- uniq

# Set 1 if the SNP appears in the sample (otherwise the matrix is 0)
ix <- unique(cbind(as.character(data$Sample), as.character(data$Key)))
m[ix] <- 1

d <- proxy::dist(m, method='Jaccard')
h = hclust(d)
plot(h)
