suppressMessages(library(stringr))
suppressMessages(library(devtools))
suppressMessages(library(VariantAnnotation))
suppressMessages(install_github("PapenfussLab/StructuralVariantAnnotation"))
suppressMessages(library(StructuralVariantAnnotation))

args <- commandArgs(trailingOnly = TRUE)
print(paste("Loading", args[[1]]))

vcf <- readVcf(args[[1]], "hg19")

normal <- rownames(colData(vcf))[[1]] # Assume the first sample is normal
tumor  <- rownames(colData(vcf))[[2]] # Assume the second sample is tumor

print(paste('Normal:', normal))
print(paste('Tumor:', tumor))

# Filter out low quality calls
vcf <- vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS"),]

# Somatic calls have no support in the normal
somatic_vcf <- vcf[geno(vcf)$QUAL[, normal] == 0,]

simpleEventType <- function(gr) {
    return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
                  ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
                         ifelse(strand(gr) == strand(partner(gr)), "INV",
                                ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
                                       "DUP")))))
}

# Extract the breakpoints by the author's package
gr <- breakpointRanges(somatic_vcf)

bed <- data.frame(chrom1=seqnames(gr),
                  start1=start(gr) - 1,
                  end1=end(gr),
                  chrom1=seqnames(partner(gr)),
                  start2=start(partner(gr)) - 1,
                  end2=end(partner(gr)),
                  name=names(gr),
                  score=gr$QUAL,
                  type=simpleEventType(gr),
                  strand1=strand(gr),
                  strand2=strand(partner(gr)))

# Just the lower of the two breakends so we don't output everything twice
bed <- bed[str_detect(bed$name, "gridss.+o"),]

print(paste("Writing to", args[[2]]))
write.table(bed, args[[2]], quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
