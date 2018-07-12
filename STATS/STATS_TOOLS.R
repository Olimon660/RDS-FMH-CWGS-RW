
parseFASTQC <- function(sample)
{
    f1 <- paste('3/FASTQC_RESULTS/F1_PAIRED_', sample, '_fastqc/fastqc_data.txt', sep='')
    f2 <- paste('3/FASTQC_RESULTS/F2_PAIRED_', sample, '_fastqc/fastqc_data.txt', sep='')
    F1_system_out <- system(paste("awk '/>>Per sequence quality scores/{f=1;next}/>>END_MODULE/{f=0}f'", f1), intern = TRUE)
    F2_system_out <- system(paste("awk '/>>Per sequence quality scores/{f=1;next}/>>END_MODULE/{f=0}f'", f2), intern = TRUE)
    F1_list <- str_split(F1_system_out[-1],"\t") # remove header
    F2_list <- str_split(F2_system_out[-1],"\t") # remove header

    df1 <- do.call(rbind.data.frame, F1_list)
    df2 <- do.call(rbind.data.frame, F2_list)
    colnames(df1) <- c('Phred', 'Count')
    colnames(df2) <- c('Phred', 'Count')
    df1$Phred <- as.numeric(as.character(df1$Phred))
    df2$Phred <- as.numeric(as.character(df2$Phred))    
    df1$Count <- as.numeric(as.character(df1$Count))
    df2$Count <- as.numeric(as.character(df2$Count))    
    
    F1_Total <- sum(df1$Count)
    F2_Total <- sum(df2$Count)
    F1_Q30_count <-sum(df1$Count[df1$Phred >=30])
    F2_Q30_count <-sum(df2$Count[df2$Phred >=30])
    F1_percent_Q30<- F1_Q30_count/F1_Total
    F2_percent_Q30<- F2_Q30_count/F2_Total
    
    data.frame(Sample=sample, F1Total=F1_Total, F2Total=F2_Total, F1Q30=F1_percent_Q30, F2Q30=F2_percent_Q30)
}

parseSAMFlag <- function(file)
{
    samtools_list                             <- scan(file, sep="\n", what=character())
    samtools_list                             <- str_trim(gsub("\\([^()]*\\)", "", samtools_list))
    passed_reads                              <- unlist(str_split(samtools_list[1],' '))
    duplicates                                <- unlist(str_split(samtools_list[2],' '))
    mapped                                    <- unlist(str_split(samtools_list[3],' '))
    paired                                    <- unlist(str_split(samtools_list[4],' '))
    read1                                     <- unlist(str_split(samtools_list[5],' '))
    read2                                     <- unlist(str_split(samtools_list[6],' '))
    properly_paired                           <- unlist(str_split(samtools_list[7],' '))
    itself_and_mate_mapped                    <- unlist(str_split(samtools_list[8],' '))
    singletons                                <- unlist(str_split(samtools_list[9],' '))
    with_mate_mapped_to_a_different_chr       <- unlist(str_split(samtools_list[10],' '))
    with_mate_mapped_to_a_different_chr_mapQ5 <- unlist(str_split(samtools_list[11],' '))
    
    data.frame(File=file,
               total_reads=passed_reads[1],
               duplicate_reads=duplicates[1], 
               mapped_reads=mapped[1], 
               all_paired_reads=paired[1], 
               count_reads1=read1[1],
               count_reads2=read2[1],
               proper_pair=properly_paired[1], 
               with_itself_and_mate_mapped=itself_and_mate_mapped[1], 
               singletons=singletons[1], 
               with_mate_mapped_to_a_different_chr= with_mate_mapped_to_a_different_chr[1],
               with_mate_mapped_to_a_different_chr_mapQ5=with_mate_mapped_to_a_different_chr_mapQ5[1])
}

parseCollectMultipleMetricsInsertSizeMetrics <- function(file)
{
    system(paste('cat', file, '| grep -v "^#" | head -n3 | tail -n2 > /tmp/A.txt', sep=' '))
    insert_size_read <- read.table(file.path('/tmp/A.txt'), sep="\t", header=TRUE)
    data.frame(File=file, InsertSize=insert_size_read)
}

parseCollectMultipleMetricsAlign <- function(file)
{
    system(paste('cat', file, '| grep -v "^#" > /tmp/A.txt', sep=' '))
    x <- read.table(file.path('/tmp/A.txt'), skip=1, sep="\t", header=TRUE)
    x <- x[x$CATEGORY == 'PAIR',]
    data.frame(File=file, AligmentRead=x)
}

parseWGS <- function(file)
{
    system(paste('cat', file, '| grep -v "^#" | head -3 > /tmp/A.txt', sep=' '))    
    data.frame(File=file, read.table(file.path('/tmp/A.txt'), skip=1, sep="\t", header=TRUE))
}

parseTrim <- function(file)
{
    out <- system(paste('cat', file, '| grep \"Input Read Pairs\"', sep=' '), intern=TRUE)
    
    if (length(out)> 0)
    { 
        trimmed_line       <- str_trim(gsub("\\([^()]*\\)", "", out))
        array              <- unlist(str_split(trimmed_line,' '))
        input_read_pairs   <- as.numeric(as.character(array[4]))
        both_surviving     <- as.numeric(as.character(array[7]))
        PCT_both_surviving <- both_surviving/input_read_pairs 
        
        forward_only_surviving     <- as.numeric(as.character(array[12]))
        PCT_forward_only_surviving <- forward_only_surviving/input_read_pairs 
        
        reverse_only_surviving     <- as.numeric(as.character(array[17]))
        PCT_reverse_only_surviving <- reverse_only_surviving/input_read_pairs
        
        dropped     <- as.numeric(as.character(array[20]))
        PCT_dropped <- dropped/input_read_pairs 
        
        x <- data.frame(File=file, InputReadPairs=input_read_pairs, BothSurviving=both_surviving, PCTBothSurviving=PCT_both_surviving, ForwardOnlySurviving=forward_only_surviving,
                        PCTForwardOnlySurviving=PCT_forward_only_surviving, ReverseOnlySurviving=reverse_only_surviving, PCTReverseOnlySurviving=PCT_reverse_only_surviving, Dropped=dropped, PCTDropped=PCT_dropped)
    }
}



