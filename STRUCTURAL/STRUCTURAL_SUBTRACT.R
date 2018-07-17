#suppressMessages(library(stringr))
#suppressMessages(library(devtools))
#suppressMessages(library(VariantAnnotation))
#suppressMessages(install_github("PapenfussLab/StructuralVariantAnnotation"))
#suppressMessages(library(StructuralVariantAnnotation))

sample <- function(x)
{
    gsub("8/ANNOTATED_SOMATIC_GRIDSS_", "", gsub(".vcf", "", x))
}

data1 <- NULL
data2 <- NULL

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
    
    samp <- sample(file)
    writeVcf(somatic_vcf, filename=paste('8/', 'SUBTRACTED_', basename(file), sep=''))
    
    simpleEventType <- function(gr) {
        return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
               ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
               ifelse(strand(gr) == strand(partner(gr)), "INV",
               ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL", "DUP")))))
    }
    
    suppressWarnings(b <- breakpointRanges(somatic_vcf))
    
    tmp1 <- data.frame(SVType  = simpleEventType(b),
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
                       Sample  = samp)
    tmp1 <- tmp1[str_detect(tmp1$Name, "gridss.+o"),] # Just the lower of the two breakends so we don't output everything twice
    
    i1 <- info(somatic_vcf[i$PARID %in% row.names(somatic_vcf),])

    
    ab <- somatic_vcf[i$PARID %in% row.names(somatic_vcf),]
    
    i$PARID
        i2 <- info(partner(somatic_vcf))
    
    tmp2 <- data.frame(Resolution = ifelse(i$IMPRECISE, "imprecise", "precise"),
                       Chr1       = seqnames(somatic_vcf),
                       Start1     = start(somatic_vcf),
                       End1       = end(somatic_vcf),
                       Chr2       = seqnames(partner(somatic_vcf)),
                       Start2     = start(partner(somatic_vcf)),
                       End2       = end(partner(somatic_vcf)),
                       Ref1       = ref(somatic_vcf),
                       Alt1       = alt(somatic_vcf),
                       Ref2       = ref(partner(somatic_vcf)),
                       Alt2       = alt(partner(somatic_vcf)),
                       RP1        = i1$RP,
                       RP2        = somatic_vcf[partner(somatic_vcf) %in% row.names(somatic_vcf),]
                       
                       )

    #data1 <- rbind(data1, ) 
}

write.table(data1, file='8/STRUCTURAL_SUBTRACT.tsv', sep='\t', quote=F, row.names=F)
