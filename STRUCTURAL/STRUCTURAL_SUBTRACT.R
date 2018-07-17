#suppressMessages(library(stringr))
#suppressMessages(library(devtools))
#suppressMessages(library(VariantAnnotation))
#suppressMessages(install_github("PapenfussLab/StructuralVariantAnnotation"))
#suppressMessages(library(StructuralVariantAnnotation))

sample <- function(x)
{
    gsub("8/ANNOTATED_SOMATIC_GRIDSS_", "", gsub(".vcf", "", x))
}

data <- NULL

for (file in Sys.glob("8/*ANNOTATED_*vcf"))
{
    print(paste("Loading", file))
    vcf <- readVcf(file, "hg19")

    normal <- rownames(colData(vcf))[[1]] # Assume the first sample is normal
    tumor  <- rownames(colData(vcf))[[2]] # Assume the second sample is tumor
    
    print(paste('Normal:', normal))
    print(paste('Tumor:', tumor))
    
    # Somatic calls have no support in the normal
    somatic_vcf <- vcf[geno(vcf)$QUAL[,normal] == 0,]
    
    simpleEventType <- function(gr) {
        return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
               ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
               ifelse(strand(gr) == strand(partner(gr)), "INV",
               ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL", "DUP")))))
    }
    
    suppressWarnings(b <- breakpointRanges(somatic_vcf))
    
    df2 <- data.frame(SVType  = simpleEventType(b),
                      Chr1    = seqnames(b),
                      Start1  = start(b) - 1,
                      End1    = end(b),
                      Chr2    = seqnames(partner(b)),
                      Start2  = start(partner(b)) - 1,
                      End2    = end(partner(b)),
                      Ref1    = b$REF,
                      Alt1    = b$ALT,
                      Name    = names(b),
                      Partner = names(partner(b)),
                      Ref2    = partner(b)$REF,
                      Alt2    = partner(b)$ALT,
                      SVLen   = b$svLen,
                      File    = file,
                      Sample  = sample(file))
    
    # Just the lower of the two breakends so we don't output everything twice
    data <- rbind(data, df2[str_detect(df2$Name, "gridss.+o"),])
}

write.table(data, file='8/STRUCTURAL_SUBTRACT.tsv', sep='\t', quote=F)
