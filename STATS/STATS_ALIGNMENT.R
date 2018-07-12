library(stringr)
source("STATS/STATS_TOOLS.R")

wgs <- data.frame()
for (file in list.files('9', pattern=glob2rx('*WGS*'))) {
    wgs <- rbind(wgs, parseWGS(paste('9/', file, sep='')))
}

samtools <- data.frame()
for (file in list.files('9', pattern=glob2rx('SAMTOOLS_*.txt'))) {
    samtools <- rbind(samtools, parseSAMFlag(paste('9/', file, sep='')))
}

fastqc <- data.frame()
for (path in list.dirs("3/FASTQC_RESULTS", full.names=FALSE, recursive=FALSE)) {
    fastqc <- rbind(fastqc, parseFASTQC(paste("3/FASTQC_RESULTS", path, 'fastqc_data.txt', sep='/')))
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
    keys <- c('9/COLLECT_WGS_METRICS_REALIGNED_RG_DEDUP_SORTED_HG19_', '.txt', '.alignment_summary_metrics', 'FLAGSTAT_REALIGNED_RG_DEDUP_SORTED_', '3/FASTQC_RESULTS/F1_PAIRED_', '3/FASTQC_RESULTS/F2_PAIRED_', '_fastqc/fastqc_data.txt', '.bam', '.insert_size_metrics', 'SAMTOOLS_FLAGSTAT_REALIGNED_RG_DEDUP_SORTED_HG19_', 'FASTQC_RESULTS/F2_PAIRED_', '.html', '2/TRIMMOMATIC_', '.log', '9/COLLECT_MULTIPLE_METRICS_REALIGNED_RG_DEDUP_SORTED_HG19_', '.bam.alignment_summary_metrics')
    for (key in keys) { x <- gsub(key, '', x) }
    x
}

wgs$Sample <- sample(wgs$File)
trim$Sample <- sample(trim$File)
fastqc$Sample <- sample(fastqc$File)
samtools$Sample <- sample(samtools$File)
collectMultipleMetricsAlign$Sample <- sample(collectMultipleMetricsAlign$File)
collectMultipleMetricsInsert$Sample <- sample(collectMultipleMetricsInsert$File)



