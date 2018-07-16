library(VariantAnnotation)

path <- '/Users/tedwong/Sources/RDS-FMH-CWGS-RW'
setwd(path)

data <- NULL
for (file in list.files('8', pattern=glob2rx('SUBTRACTED_*.vcf')))
{
    bed <- read.table(paste('8', file, sep='/'))
    colnames(bed) <- c('Chr1', 'Start1', 'End1', 'Chr2', 'Start2', 'End2', 'Name', 'Qual', 'Strand1', 'Strand2')
    data <- rbind(data, data.frame(File=file,
                                   Chr1=bed$Chr1,
                                   Start1=bed$Start1,
                                   End1=bed$End1,
                                   Chr2=bed$Chr2,
                                   Start2=bed$Start2,
                                   End2=bed$End2,
                                   Strand1=bed$Strand1,
                                   Strand2=bed$Strand2))
}
