suppressMessages(library(stringr))
suppressMessages(library(devtools))
suppressMessages(library(VariantAnnotation))
suppressMessages(install_github("PapenfussLab/StructuralVariantAnnotation"))
suppressMessages(library(StructuralVariantAnnotation))

sample <- function(x) {
    gsub("8/ANNOTATED_SOMATIC_GRIDSS_", "", gsub(".vcf", "", x))
}

gData <- NULL
sData <- NULL

for (file in Sys.glob("8/ANNOTATED_*vcf"))
{
    print(paste("Loading", file))
    vcf <- readVcf(file, "hg19")

    normal <- rownames(colData(vcf))[[1]] # Assume the first sample is normal
    tumor  <- rownames(colData(vcf))[[2]] # Assume the second sample is tumor
    
    print(paste('Normal:', normal))
    print(paste('Tumor:', tumor))
    
    gv_ <- vcf[geno(vcf)$QUAL[,normal] != 0,]

    samp <- sample(file)
    
    simpleEventType <- function(gr) {
        return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
               ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
               ifelse(strand(gr) == strand(partner(gr)), "INV",
               ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL", "DUP")))))
    }
    
    process <- function(v) {
        # Remove unpaired breakpoints (e.g. removed by normal)
        suppressWarnings(b <- breakpointRanges(v))

        v <- v[row.names(v) %in% names(b) ,]
        writeVcf(v, filename=paste('8/', 'SUBTRACTED_', basename(file), sep=''))
        
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
        
        i1 <- info(v)
        i1 <- data.frame(Name    = row.names(v),
                         Resolu  = ifelse(i1$IMPRECISE, 'Imprecise', 'Precise'),
                         Partner = i1$PARID,
                         Filter  = filt(v),
                         RP      = i1$RP,
                         RPQ     = i1$RPQ,
                         SR      = i1$SR,
                         SRQ     = i1$SRQ,
                         REFC    = i1$REF,
                         HOMLEN  = as.numeric(i1$HOMLEN),
                         HOMSEQ  = as.character(i1$HOMSEQ),
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
                           HOMLEN1  = i2$HOMLEN.x,
                           HOMSEQ1  = i2$HOMSEQ.x,
                           HOMLEN2  = i2$HOMLEN.y,
                           HOMSEQ2  = i2$HOMSEQ.y,
                           REFPAIR1 = i2$REFPAIR.x,
                           REFPAIR2 = i2$REFPAIR.y)
        tmp3 <- merge(tmp1, tmp2, by.x="Name", by.y="Name")
        tmp3
    }
    
    gData <- rbind(gData, process(vcf[geno(vcf)$QUAL[,normal] != 0,]))
    sData <- rbind(sData, process(vcf[geno(vcf)$QUAL[,normal] == 0,])) # Overwrite the VCF file from last step
}

gData$Partner <- gData$Partner.x
sData$Partner <- sData$Partner.x
gData <- gData[,!(names(gData) %in% c("Partner.x", "Partner.y"))]
sData <- sData[,!(names(sData) %in% c("Partner.x", "Partner.y"))]

gData$Type <- "Germline"
sData$Type <- "Somatic"

write.table(rbind(gData, sData), file='8/STRUCTURAL_SUBTRACT.tsv', sep='\t', quote=F, row.names=F)
