# example sciprt performing simple filtering of variants to somatic calls
#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
library(stringr)
library(VariantAnnotation)
library(devtools)
install_github("PapenfussLab/StructuralVariantAnnotation")
library(StructuralVariantAnnotation)

vcf <- readVcf("/Users/twong/Sources/RDS-FMH-CWGS-RW/8/SOMATIC_GRIDSS_IIICF-E6_A1.vcf", "hg19")

# Filter out low quality calls
vcf <- vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS"),]

# Somatic calls have no support in the normal
somatic_vcf <- vcf[geno(vcf)$QUAL[,"SAMPLE_IIICF"] == 0,]

loh_vcf <- vcf[geno(vcf)$QUAL[,"tumour.bam"] == 0,]

# Output BEDPE for use by circos
gr <- breakpointRanges(somatic_vcf)
bedpe <- data.frame(
    chrom1=seqnames(gr),
    start1=start(gr) - 1,
    end1=end(gr),
    chrom1=seqnames(partner(gr)),
    start2=start(partner(gr)) - 1,
    end2=end(partner(gr)),
    name=names(gr),
    score=gr$QUAL,
    strand1=strand(gr),
    strand2=strand(partner(gr))
)

# Just the lower of the two breakends so we don't output everything twice
bedpe <- bedpe[str_detect(bedpe$name, "gridss.+o"),]
write.table(bedpe, "somatic.gridss.hq.somatic.bedpe", quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)