library(dplyr)
library(stringr)
source("STATS/STATS_TOOLS.R")

x <- list.dirs("3/FASTQC_RESULTS", full.names=FALSE, recursive=FALSE)
x <- gsub('F1_', '', x)
x <- gsub('F2_', '', x)
x <- gsub('_fastqc', '', x)
x <- gsub('PAIRED_', '', x)
x <- unique(x)

fastqc <- data.frame()
for (sample in x) {
    fastqc <- rbind(fastqc, parseFASTQC(sample))
}

wgs <- data.frame()
for (file in list.files('9', pattern=glob2rx('*WGS*'))) {
    wgs <- rbind(wgs, parseWGS(paste('9/', file, sep='')))
}

samtools <- data.frame()
for (file in list.files('9', pattern=glob2rx('SAMTOOLS_*.txt'))) {
    samtools <- rbind(samtools, parseSAMFlag(paste('9/', file, sep='')))
}

collectMultipleMetricsInsert <- data.frame()
for (file in gsub('.txt', '', list.files('9', pattern=glob2rx('*insert_size_metrics')))) {
    collectMultipleMetricsInsert <- rbind(collectMultipleMetricsInsert, parseCollectMultipleMetricsInsertSizeMetrics(paste('9/', file, sep='')))
}

collectMultipleMetricsAlign <- data.frame()
for (file in gsub('.txt', '', list.files('9', pattern=glob2rx('*alignment_summary_metrics')))) {
    collectMultipleMetricsAlign <- rbind(collectMultipleMetricsAlign, parseCollectMultipleMetricsAlign(paste('9/', file, sep='')))
}

trim <- data.frame()
for (file in list.files('2', pattern=glob2rx('*.log'))) {
    trim <- rbind(trim, parseTrim(paste('2/', file, sep='')))
}

sample <- function(x)
{
    keys <- c('9/SAMTOOLS_FLAGSTAT_REALIGNED_RG_DEDUP_SORTED_HG19_', '9/SAMTOOLS_HG19_', '9/COLLECT_WGS_METRICS_REALIGNED_RG_DEDUP_SORTED_HG19_', '.txt', '.alignment_summary_metrics', 'FLAGSTAT_REALIGNED_RG_DEDUP_SORTED_', '3/FASTQC_RESULTS/F1_PAIRED_', '3/FASTQC_RESULTS/F2_PAIRED_', '_fastqc/fastqc_data.txt', '.bam', '.insert_size_metrics', 'SAMTOOLS_FLAGSTAT_REALIGNED_RG_DEDUP_SORTED_HG19_', 'FASTQC_RESULTS/F2_PAIRED_', '.html', '2/TRIMMOMATIC_', '.log', '9/COLLECT_MULTIPLE_METRICS_REALIGNED_RG_DEDUP_SORTED_HG19_', '.bam.alignment_summary_metrics')
    for (key in keys) { x <- gsub(key, '', x) }
    x
}

wgs$Sample <- sample(wgs$File)
trim$Sample <- sample(trim$File)
samtools$Sample <- sample(samtools$File)
collectMultipleMetricsAlign$Sample <- sample(collectMultipleMetricsAlign$File)
collectMultipleMetricsInsert$Sample <- sample(collectMultipleMetricsInsert$File)

# No trimming
specials1 <- c('H06L4ALXX_2_JFCF6_P-pLKO_5_Human__Reddel_Lab', 'JFCF6_P12_H06L4ALXX_5', 'JFCF6_T_1J_1-3C_H0AF0ALXX_3', 'JFCF6_T_1J_6B_H0AF0ALXX_4')

# No FASTQC
specials2 <- c('IIICF', 'JFCF_61Q')

wgs  <- wgs[!(wgs$Sample %in% specials1) & !(wgs$Sample %in% specials2),]
trim <- trim[!(trim$Sample %in% specials1) & !(trim$Sample %in% specials2),]
fastqc <- fastqc[!(fastqc$Sample %in% specials1) & !(fastqc$Sample %in% specials2),]
samtools <- samtools[!(samtools$Sample %in% specials1) & !(samtools$Sample %in% specials2),]
collectMultipleMetricsInsert <- collectMultipleMetricsInsert[!(collectMultipleMetricsInsert$Sample %in% specials1) & !(collectMultipleMetricsInsert$Sample %in% specials2),]
collectMultipleMetricsAlign <- collectMultipleMetricsAlign[!(collectMultipleMetricsAlign$Sample %in% specials1) & !(collectMultipleMetricsAlign$Sample %in% specials2),]

wgs <- dplyr::select(wgs, -c(File))
trim <- dplyr::select(trim, -c(File))
samtools <- dplyr::select(samtools, -c(File))
collectMultipleMetricsAlign <- dplyr::select(collectMultipleMetricsAlign, -c(File))
collectMultipleMetricsInsert <- dplyr::select(collectMultipleMetricsInsert, -c(File))

sort <- function(x) { x[with(x, order(Sample)), ] }

wgs <- sort(wgs)
trim <- sort(trim)
samtools <- sort(samtools)
collectMultipleMetricsAlign <- sort(collectMultipleMetricsAlign)
collectMultipleMetricsInsert <- sort(collectMultipleMetricsInsert)

colnames(samtools) <- c('Raw Total Individual Reads')

data <- data.frame(SampleName=samtools$Sample, TotalReads=samtools$total_reads)




