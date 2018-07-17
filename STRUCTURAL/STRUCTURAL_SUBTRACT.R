suppressMessages(library(stringr))
suppressMessages(library(devtools))
suppressMessages(library(VariantAnnotation))
suppressMessages(install_github("PapenfussLab/StructuralVariantAnnotation"))
suppressMessages(library(StructuralVariantAnnotation))

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
    
    samp <- sample(file)
    
    simpleEventType <- function(gr) {
        return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
               ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
               ifelse(strand(gr) == strand(partner(gr)), "INV",
               ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL", "DUP")))))
    }
    
    # Remove unpaired breakpoints (e.g. removed by normal)
    suppressWarnings(b <- breakpointRanges(somatic_vcf))
    
    somatic_vcf <- somatic_vcf[row.names(somatic_vcf) %in% names(b) ,]
    writeVcf(somatic_vcf, filename=paste('8/', 'SUBTRACTED_', basename(file), sep=''))

    # Information from the SV package
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
    
    i1 <- info(somatic_vcf)
    
    # Information from the VCF
    i1 <- data.frame(Name    = row.names(somatic_vcf),
                     Resolu  = ifelse(i1$IMPRECISE, 'Imprecise', 'Precise'),
                     Partner = i1$PARID,
                     Filter  = filt(somatic_vcf),
                     RP      = i1$RP,
                     RPQ     = i1$RPQ,
                     SR      = i1$SR,
                     SRQ     = i1$SRQ,
                     REFC    = i1$REF,
                     REFPAIR = i1$REFPAIR)
    i1_ <- i1
    i1_ <- i1_[,!(names(i1_) %in% c("Name"))]
    i2   <- merge(x=i1, y=i1_, by.x="Name", by.y="Partner")
    tmp2 <- data.frame(Name     = i2$Name,
                       Resolu   = i2$Resolu.x,
                       Partner  = i2$Partner,
                       Filter   = i2$Filter.x,
                       RP1      = i2$RP.x,
                       RP2      = i2$RP.y,
                       RPQ1     = i2$RPQ.x,
                       RPQ2     = i2$RPQ.y,
                       SR1      = i2$SR.x,
                       SR2      = i2$SR.y,
                       SRQ1     = i2$SRQ.x,
                       SRQ2     = i2$SRQ.y,
                       REFC1    = i2$REFC.x,
                       REFC2    = i2$REFC.y,
                       REFPAIR1 = i2$REFPAIR.x,
                       REFPAIR2 = i2$REFPAIR.y)
    
    tmp3 <- merge(tmp1, tmp2, by.x="Name", by.y="Name")
    data <- rbind(data, tmp3)
}

data$Partner <- data$Partner.x
data <- data[,!(names(data) %in% c("Partner.x", "Partner.y"))]
write.table(data, file='8/STRUCTURAL_SUBTRACT.tsv', sep='\t', quote=F, row.names=F)
