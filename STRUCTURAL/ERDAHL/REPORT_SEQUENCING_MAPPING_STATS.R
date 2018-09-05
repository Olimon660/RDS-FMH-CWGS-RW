
temp1<-list.files('/mnt/TempSAN/eteber/Macrogen_Sequenced', '*.fq')
temp1<-gsub('F1_PAIRED_','',temp1)
temp1<-gsub('F2_PAIRED_','',temp1)
temp1<-gsub('F1_UNPAIRED_','',temp1)
temp1<-gsub('F2_UNPAIRED_','',temp1)
temp1<-gsub('.fq','',temp1)
temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

script<-c()

command <-'/usr/bin/bioinformatics/FastQC/fastqc'
output_dir<-'--outdir=/mnt/TempSAN/eteber/Macrogen_Sequenced/cleaned_FASTQ_QC'

for (sampleID in sample_list)
{
  
  #/usr/bin/bioinformatics/FastQC/fastqc  /mnt/san2/lab1/Bioinfo_Projects/Garvan_WGS_OCT_2014/2_clean_FASTQ/F1_PAIRED_H0A1UALXX_2_CHLA_291_C3.fq  --outdir=/mnt/san2/lab1/Bioinfo_Projects/Garvan_WGS_OCT_2014/3_cleaned_FASTQ_QC 
  #/usr/bin/bioinformatics/FastQC/fastqc  /mnt/san2/lab1/Bioinfo_Projects/Garvan_WGS_OCT_2014/2_clean_FASTQ/F2_PAIRED_H0A1UALXX_2_CHLA_291_C3.fq  --outdir=/mnt/san2/lab1/Bioinfo_Projects/Garvan_WGS_OCT_2014/3_cleaned_FASTQ_QC 
  
  print (sampleID)
  
  f1_paired_filename<-paste('F1_PAIRED_', sampleID, '.fq', sep='')
  f1_paired_input   <-file.path('/mnt/TempSAN/eteber/Macrogen_Sequenced',f1_paired_filename)
  
  f2_paired_filename<-paste('F2_PAIRED_', sampleID, '.fq', sep='')
  f2_paired_input   <-file.path('/mnt/TempSAN/eteber/Macrogen_Sequenced',f2_paired_filename)
  
  
  script<- c(script, (paste(command,f1_paired_input, output_dir,  sep=' ')))
  script<- c(script, (paste(command,f2_paired_input, output_dir,  sep=' ')))
  
  
}

write(script, '/mnt/san2/bioinformatics/eteber/cancer_group/rdagg/CLEANED_FASTQC_SCRIPT.txt')


parallel -j 6 < /mnt/san2/bioinformatics/eteber/cancer_group/rdagg/CLEANED_FASTQC_SCRIPT.txt




####################################################################################################################################



script<-c()
temp1<-list.files('/mnt/san2/lab1/Bioinfo_Projects/mlee_marked_dups', pattern=glob2rx('DEDUP_SORTED_*.bam'))
temp1<-gsub('.bam','',temp1)
temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

command   <-'java  -jar /mnt/san2/bioinformatics/Applications/picard-tools-1.120/CollectMultipleMetrics.jar'
input_dir <-'/mnt/san2/lab1/Bioinfo_Projects/mlee_marked_dups'
output_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'
additional_parameters<-'R=/mnt/san2/bioinformatics/eteber/hg19/ucsc.hg19.fasta TMP_DIR=/mnt/san2/tmp/eteber'
programs<-'PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle'

library(stringr) 

for (sampleID in sample_list)
{
  
  string_list<-unlist(str_split(sampleID,'_'))
  sample_name<-paste(string_list[1:length(string_list)],collapse = '_')
  out_filename<-paste('CollectMultipleMetrics_', sample_name,'.txt',sep='')
  print (sampleID)
  
  input_string  <-paste('INPUT=', file.path(input_dir, sampleID), '.bam',sep='')
  output_string <-paste('OUTPUT=', file.path(output_dir, out_filename),sep='')
  error_string  <-paste('2> ', file.path(output_dir, paste(sample_name, '_CollectMultipleMetrics.err', sep='')), sep='')
  
  script<- c(script, paste(command, input_string, output_string, programs, additional_parameters, sep=' '))
  
}

write(script, '/mnt/san2/bioinformatics/eteber/cancer_group/rdagg/CollectMultipleMetrics.txt')


parallel -j 10 < /mnt/san2/bioinformatics/eteber/cancer_group/rdagg/CollectMultipleMetrics.txt


TMP=/mnt/san2/temp
TMPDIR=/mnt/san2/temp
export TMP TMPDIR


#/mnt/san2/bioinformatics/Applications/picard-tools-1.120/CollectWgsMetrics.jar

java  -jar /mnt/san2/bioinformatics/Applications/picard-tools-1.120/CollectWgsMetrics.jar 
INPUT=/mnt/san2/lab1/Bioinfo_Projects/mlee_marked_dups/DEDUP_SORTED_H06L2ALXX_1_LAN_6_MUT_P375.bam 
OUTPUT=/mnt/san2/lab1/Bioinfo_Projects/mlee_marked_dups/WGS_METRICS_DEDUP_SORTED_H06L2ALXX_1_LAN_6_MUT_P375.tsv 
REFERENCE_SEQUENCE=/mnt/san2/bioinformatics/eteber/hg19/ucsc.hg19.fasta TMP_DIR=/mnt/san2/tmp/eteber

#remove histogram info from the output
#cat wgs_metric.tsv  | grep -v "^#" | head -3 > wgs_metric_no_hist.tsv

#java -jar /usr/bin/bioinformatics/picard-tools-1.119/CollectAlignmentSummaryMetrics.jar 
#INPUT=/mnt/san2/lab1/Bioinfo_Projects/Garvan_WGS_OCT_2014_Bioinformatics/4_Aligned/H06L2ALXX_5_CHLA_291_C3/NODUP_SORTED_H06L2ALXX_5_CHLA_291_C3.bam   
#OUTPUT=/mnt/san2/lab1/Bioinfo_Projects/Garvan_WGS_OCT_2014_Bioinformatics/4_Aligned/H06L2ALXX_5_CHLA_291_C3/ALIGN_METRICS_NODUP_SORTED_H06L2ALXX_5_CHLA_291_C3.tsv   
#TMP_DIR=/mnt/san2/tmp/eteber


#samtools flagstat /mnt/san2/lab1/Bioinfo_Projects/Garvan_WGS_OCT_2014_Bioinformatics/4_Aligned/H06L2ALXX_5_CHLA_291_C3/NODUP_SORTED_H06L2ALXX_5_CHLA_291_C3.bam > /mnt/san2/lab1/Bioinfo_Projects/Garvan_WGS_OCT_2014_Bioinformatics/4_Aligned/H06L2ALXX_5_CHLA_291_C3/SAMTOOLS_NODUP_SORTED_H06L2ALXX_5_CHLA_291_C3.txt

script<-c()
temp1<-list.files('/mnt/san2/lab1/Bioinfo_Projects/mlee_marked_dups', pattern=glob2rx('DEDUP_SORTED_*.bam'))
temp1<-gsub('.bam','',temp1)
temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

command   <-'java -Djava.io.tmpdir=/mnt/san2/temp/eteber -jar /mnt/san2/bioinformatics/Applications/picard-tools-1.120/CollectWgsMetrics.jar'
input_dir <-'/mnt/san2/lab1/Bioinfo_Projects/mlee_marked_dups'
output_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'

library(stringr) 

for (sampleID in sample_list)
{
  
  string_list<-unlist(str_split(sampleID,'_'))
  sample_name<-paste(string_list[1:length(string_list)],collapse = '_')
  out_filename<-paste('WGS_METRICS_', sample_name,'.tsv',sep='')
  print (sampleID)
  
  input_string  <-paste('INPUT=', file.path(input_dir, sampleID), '.bam',sep='')
  output_command<-paste('OUTPUT=',file.path(output_dir, out_filename), sep='')
  tail_command<-'REFERENCE_SEQUENCE=/mnt/san2/bioinformatics/eteber/hg19/ucsc.hg19.fasta TMP_DIR=/mnt/san2/tmp/eteber'
  
  script<- c(script, paste(command, input_string, output_command, tail_command, sep=' '))
  
}

write(script, '/mnt/san2/bioinformatics/eteber/cancer_group/rdagg/CollectWgsMetrics.txt')


parallel -j 22 < /mnt/san2/bioinformatics/eteber/cancer_group/rdagg/CollectWgsMetrics.txt


#remove histogram info from the output
#cat wgs_metric.tsv  | grep -v "^#" | head -3 > wgs_metric_no_hist.tsv

script<-c()
temp1<-list.files('/mnt/san2/lab1/Bioinfo_Projects/mlee_marked_dups', pattern=glob2rx('DEDUP_SORTED_*.bam'))
temp1<-gsub('.bam','',temp1)
temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

command   <-'samtools flagstat'
input_dir <-'/mnt/san2/lab1/Bioinfo_Projects/mlee_marked_dups'
output_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'

library(stringr) 

for (sampleID in sample_list)
{
  
  string_list<-unlist(str_split(sampleID,'_'))
  sample_name<-paste(string_list[1:length(string_list)],collapse = '_')
  out_filename<-paste('SAMTOOLS_', sample_name,'.txt',sep='')
  print (sampleID)
  
  input_string  <-paste(file.path(input_dir, sampleID), '.bam',sep='')
  output_command<-paste('> ',file.path(output_dir, out_filename), sep='')
  
  script<- c(script, paste(command, input_string, output_command, sep=' '))
  
}

write(script, '/mnt/san2/bioinformatics/eteber/cancer_group/rdagg/SAMTOOLSMetrics.txt')


parallel -j 22 < /mnt/san2/bioinformatics/eteber/cancer_group/rdagg/SAMTOOLSMetrics.txt


#######################################################
#
####START HERE FOR REPORTS:
#  

# read in details


#for file in SAMTOOLS_DEDUP_SORTED_SORTED_*; do cp "$file" "${file/SAMTOOLS_DEDUP_SORTED_SORTED_/SAMTOOLS_DEDUP_SORTED_}";done


input_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'
script<-c()
temp1<-list.files('/mnt/TempSAN/eteber/Sequencing_Mapping_Stats', pattern=glob2rx('SAMTOOLS_*.txt'))
temp1<-gsub('.txt','',temp1)
temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

library(stringr) 
db_samtools <- data.frame()

for (sampleID in sample_list)
{
   
  sample_file<-paste(sampleID,'.txt',sep='')
  
  sample_name<-gsub('SAMTOOLS_','',sampleID)
  
  samtools_list                             <- scan(file.path(input_dir,sample_file), sep="\n", what=character())
  samtools_list                             <- str_trim(gsub("\\([^()]*\\)", "", samtools_list))
 
  passed_reads                              <- unlist(str_split(samtools_list[1],' '))   #763868851 + 0 in total (QC-passed reads + QC-failed reads)
  
  duplicates                                <- unlist(str_split(samtools_list[2],' '))   #169544006 + 0 duplicates
  
  mapped                                    <- unlist(str_split(samtools_list[3],' '))   #763134788 + 0 mapped (99.90%:-nan%)
  
  paired                                    <- unlist(str_split(samtools_list[4],' '))   #763868851 + 0 paired in sequencing 

  read1                                     <- unlist(str_split(samtools_list[5],' '))   #382019437 + 0 read1

  read2                                     <- unlist(str_split(samtools_list[6],' '))   #381849414 + 0 read2
  
  properly_paired                           <- unlist(str_split(samtools_list[7],' '))   #745766626 + 0 properly paired (97.63%:-nan%)

  itself_and_mate_mapped                    <- unlist(str_split(samtools_list[8],' '))   #762911434 + 0 with itself and mate mapped

  singletons                                <- unlist(str_split(samtools_list[9],' '))   #   223354 + 0 singletons (0.03%:-nan%)

  with_mate_mapped_to_a_different_chr       <- unlist(str_split(samtools_list[10],' '))  # 10154353 + 0 with mate mapped to a different chr

  with_mate_mapped_to_a_different_chr_mapQ5 <- unlist(str_split(samtools_list[11],' '))  #  3303463 + 0 with mate mapped to a different chr (mapQ>=5)

  my_record<-data.frame(sample_name=sample_name, total_reads=passed_reads[1],
                        duplicate_reads=duplicates[1], 
                        mapped_reads=mapped[1], 
                        all_paired_reads= paired[1], 
                        count_reads1 = read1[1],
                        count_reads2 = read2[1],
                        proper_pair=properly_paired[1], 
                        with_itself_and_mate_mapped=itself_and_mate_mapped[1], 
                        singletons=singletons[1], 
                        with_mate_mapped_to_a_different_chr= with_mate_mapped_to_a_different_chr[1],
                        with_mate_mapped_to_a_different_chr_mapQ5=with_mate_mapped_to_a_different_chr_mapQ5[1])

  
  #n_reads are the total number of reads
  #n_pair_all : the read is paired in sequencing, no matter whether it is mapped in a pair
  #n_pair_good : the read is mapped in a proper pair
  #n_read1 : count read1
  #n_read2 : count read2
  #n_sgltn : the read itself is unmapped the mate is mapped
  #n_pair_map: the read itself is mapped the mate is unmapped
  #n_diffchr: number of reads with a mate mapped on a different chromosome
  #n_diffhigh: number of reads with a mate on a different chromosome having a quality greater than 5
  
  db_samtools<-rbind(db_samtools,my_record)
  
}




###############################################################################


input_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'
script<-c()
temp1<-list.files('/mnt/TempSAN/eteber/Sequencing_Mapping_Stats', pattern=glob2rx('WGS_*.tsv'))
temp1<-gsub('.tsv','',temp1)
temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

library(stringr) 


for (sampleID in sample_list)
{
  
  sample_file     <-paste(sampleID,'.tsv',sep='')
  out_file        <-paste('SMALL_','',sample_file, sep='')
  sh_string<-paste('cat', file.path(input_dir,sample_file), '| grep -v "^#" | head -3 >', file.path(input_dir,out_file), sep=' ') 
  system(sh_string)
}

##############################

input_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'
temp1<-list.files('/mnt/TempSAN/eteber/Sequencing_Mapping_Stats', pattern=glob2rx('SMALL_WGS_*.tsv'))
temp1<-gsub('.tsv','',temp1)
temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

library(stringr) 

db_picard_wgstools <- data.frame()

for (sampleID in sample_list)
{
  sample_file       <- paste(sampleID,'.tsv',sep='')
  WGS_read          <- read.table(file.path(input_dir,sample_file), skip=1, sep="\t", header=TRUE)
  sample_name       <- gsub('SMALL_WGS_METRICS_DEDUP_SORTED__','DEDUP_SORTED_',sampleID)
  wgs_record        <- data.frame(sample_name, WGS_read)
  db_picard_wgstools<-rbind(db_picard_wgstools,wgs_record)
}

library(sqldf)




################################################################################################

#read FASTQC details

#/usr/bin/bioinformatics/FastQC/fastqc /mnt/san2/temp/eteber/F1_PAIRED_TEST.fq --outdir=/mnt/san2/temp/eteber

#cat /mnt/TempSAN/eteber/3_cleaned_FASTQ_QC/F1_PAIRED_H0A1UALXX_2_CHLA_291_C3.fq_fastqc/fastqc_data.txt | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > F1_PAIRED_H0A1UALXX_2_CHLA_291_C3_LENGTHS.txt

#>>Per base sequence quality

#>>END_MODULE




input_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'

temp1<-list.dirs('/mnt/TempSAN/eteber/Sequencing_Mapping_Stats', full.names=FALSE, recursive = FALSE)
temp1<-gsub('F1_PAIRED_','',temp1)
temp1<-gsub('F2_PAIRED_','',temp1)
temp1<-gsub('.fq_fastqc','',temp1)

temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

library(stringr) 

db_fastqc <- data.frame()

for (sampleID in sample_list)
{

  F1_sample_dir        <- paste('F1_PAIRED_',sampleID,'.fq_fastqc',sep='')
  F2_sample_dir        <- paste('F2_PAIRED_',sampleID,'.fq_fastqc',sep='')

  #awk '/>>Per sequence quality scores/{f=1;next}/>>END_MODULE/{f=0}f' /mnt/TempSAN/eteber/Sequencing_Mapping_Stats/F1_PAIRED_H0A25ALXX_7_LAN_6_wt_P259.fq_fastqc/fastqc_data.txt
  
  F1_sh_string<-paste("awk '/>>Per sequence quality scores/{f=1;next}/>>END_MODULE/{f=0}f'", file.path(input_dir, F1_sample_dir, "fastqc_data.txt"))
  F1_system_out<- system(F1_sh_string,intern = TRUE)
  F1_list<-str_split(F1_system_out[-1],"\t") # remove header
  F1_db<-do.call(rbind.data.frame, F1_list)  
  colnames(F1_db)<-c('Phred', 'Count')
  F1_db$Phred    <-as.numeric(as.character(F1_db$Phred))  
  F1_db$Count<-as.numeric(as.character(F1_db$Count))  
  F1_Q30_count  <-sum(F1_db$Count[F1_db$Phred >=30])
  F1_total_reads<-sum(F1_db$Count)
  F1_percent_Q30<- F1_Q30_count/F1_total_reads  
  
  F2_sh_string<-paste("awk '/>>Per sequence quality scores/{f=1;next}/>>END_MODULE/{f=0}f'", file.path(input_dir, F2_sample_dir, "fastqc_data.txt"))
  F2_system_out<- system(F2_sh_string,intern = TRUE)
  F2_list<-str_split(F2_system_out[-1],"\t") # remove header
  F2_db<-do.call(rbind.data.frame, F2_list)  
  colnames(F2_db)<-c('Phred', 'Count')
  F2_db$Phred    <-as.numeric(as.character(F2_db$Phred))  
  F2_db$Count    <-as.numeric(as.character(F2_db$Count))  
  F2_Q30_count  <-sum(F2_db$Count[F2_db$Phred >=30])
  F2_total_reads<-sum(F2_db$Count)
  F2_percent_Q30<- F2_Q30_count/F2_total_reads  
  
  key<-paste('DEDUP_SORTED_', sampleID, sep='')
  
  my_rec<-data.frame(key, F1_total_reads, F2_total_reads, F1_percent_Q30, F2_percent_Q30)  
  db_fastqc<-rbind(db_fastqc, my_rec)
}  


#CollectMultipleMetrics_DEDUP_SORTED_H06L2ALXX_1_LAN_6_MUT_P375.txt.insert_size_metrics

input_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'
script<-c()
temp1<-list.files('/mnt/TempSAN/eteber/Sequencing_Mapping_Stats', pattern=glob2rx('CollectMultipleMetrics_DEDUP_SORTED_*.txt.insert_size_metrics'))
temp1<-gsub('.txt.insert_size_metrics','',temp1)
temp1<-gsub('CollectMultipleMetrics_DEDUP_SORTED_','',temp1)

temp1<-unique(sort(temp1))
length(temp1)
sample_list<-temp1

library(stringr) 


for (sampleID in sample_list)
{
  
  sample_file     <-paste('CollectMultipleMetrics_DEDUP_SORTED_', sampleID,'.txt.insert_size_metrics',sep='')
  out_file        <-paste('SMALL_','',sample_file, sep='')
  sh_string<-paste('cat', file.path(input_dir,sample_file), '| grep -v "^#" | head -3 >', file.path(input_dir,out_file), sep=' ') 
  system(sh_string)
}



input_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'
temp1<-list.files('/mnt/TempSAN/eteber/Sequencing_Mapping_Stats', pattern=glob2rx('SMALL_CollectMultipleMetrics_DEDUP_SORTED_*.txt.insert_size_metrics'))
temp1<-gsub('.txt.insert_size_metrics','',temp1)
temp1<-gsub('SMALL_CollectMultipleMetrics_DEDUP_SORTED_','',temp1)
temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

library(stringr) 

db_picard_insert_size <- data.frame()

for (sampleID in sample_list)
{
  sample_file               <-paste('SMALL_CollectMultipleMetrics_DEDUP_SORTED_', sampleID,'.txt.insert_size_metrics',sep='')
  insert_size_read          <- read.table(file.path(input_dir,sample_file), skip=1, sep="\t", header=TRUE)
  sample_name               <- paste('DEDUP_SORTED_',sampleID, sep='')
  insert_size_record        <- data.frame(sample_name, insert_size_read)
  db_picard_insert_size<-rbind(db_picard_insert_size,insert_size_record)
}


#CollectMultipleMetrics_DEDUP_SORTED_SK_N_FI_H0FAV_L007.txt.alignment_summary_metrics

input_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'
script<-c()
temp1<-list.files('/mnt/TempSAN/eteber/Sequencing_Mapping_Stats', pattern=glob2rx('CollectMultipleMetrics_DEDUP_SORTED_*.txt.alignment_summary_metrics'))
temp1<-gsub('.txt.alignment_summary_metrics','',temp1)
temp1<-gsub('CollectMultipleMetrics_DEDUP_SORTED_','',temp1)

temp1<-unique(sort(temp1))
length(temp1)
sample_list<-temp1

library(stringr) 


for (sampleID in sample_list)
{
  
  sample_file     <-paste('CollectMultipleMetrics_DEDUP_SORTED_', sampleID,'.txt.alignment_summary_metrics',sep='')
  out_file        <-paste('SMALL_','',sample_file, sep='')
  sh_string<-paste('cat', file.path(input_dir,sample_file), '| grep -v "^#" >', file.path(input_dir,out_file), sep=' ') 
  system(sh_string)
}



input_dir<-'/mnt/TempSAN/eteber/Sequencing_Mapping_Stats'
temp1<-list.files('/mnt/TempSAN/eteber/Sequencing_Mapping_Stats', pattern=glob2rx('SMALL_CollectMultipleMetrics_DEDUP_SORTED_*.txt.alignment_summary_metrics'))
temp1<-gsub('.txt.alignment_summary_metrics','',temp1)
temp1<-gsub('SMALL_CollectMultipleMetrics_DEDUP_SORTED_','',temp1)
temp1<-unique(sort(temp1))
length(temp1)

sample_list<-temp1

library(stringr) 

db_picard_alignment <- data.frame()

for (sampleID in sample_list)
{
  sample_file               <-paste('SMALL_CollectMultipleMetrics_DEDUP_SORTED_', sampleID,'.txt.alignment_summary_metrics',sep='')
  alignment_read          <- read.table(file.path(input_dir,sample_file), skip=1, sep="\t", header=TRUE)
  sample_name               <- paste('DEDUP_SORTED_',sampleID, sep='')
  alignment_record        <- data.frame(sample_name, alignment_read)
  db_picard_alignment<-rbind(db_picard_alignment,alignment_record[alignment_record$CATEGORY=='PAIR',])
}

keys<-c('DEDUP_SORTED_CHLA90_H0FAV_L008',
        'DEDUP_SORTED_COG16_H0FAN_L004',
        'DEDUP_SORTED_COG21_H0FAN_L005',
        'DEDUP_SORTED_COG25_H0FAN_L006',
        'DEDUP_SORTED_COG28_H0FAN_L007',
        'DEDUP_SORTED_COG29_H0FAN_L008',
        'DEDUP_SORTED_COG7_H0FAN_L003',
        'DEDUP_SORTED_H06L2ALXX_1_LAN_6_MUT_P375',
        'DEDUP_SORTED_H06L2ALXX_2_LAN_6_pre_MUT_P81',
        'DEDUP_SORTED_H06L2ALXX_3_LAN_6_wt_P75',
        'DEDUP_SORTED_H06L2ALXX_4_LAN_6_wt_P259',
        'DEDUP_SORTED_H06L2ALXX_5_CHLA_291_C3',
        'DEDUP_SORTED_H06L2ALXX_6_CHLA_291_C4',
        'DEDUP_SORTED_H06L2ALXX_7_CHW_11',
        'DEDUP_SORTED_H06L2ALXX_8_COG_4',
        'DEDUP_SORTED_H06L4ALXX_1_COG_32',
        'DEDUP_SORTED_H09CMALXX_2_LAN_6_MUT_P375',
        'DEDUP_SORTED_H09CMALXX_3_LAN_6_pre_MUT_P81',
        'DEDUP_SORTED_H09CMALXX_4_LAN_6_wt_P75',
        'DEDUP_SORTED_H0A1UALXX_2_CHLA_291_C3',
        'DEDUP_SORTED_H0A1UALXX_3_CHLA_291_C4',
        'DEDUP_SORTED_H0A1UALXX_4_CHW_11',
        'DEDUP_SORTED_H0A1UALXX_5_COG_4',
        'DEDUP_SORTED_H0A1UALXX_6_COG_32',
        'DEDUP_SORTED_H0A25ALXX_7_LAN_6_wt_P259',
        'DEDUP_SORTED_SH_SY5Y_H0FAN_L001',
        'DEDUP_SORTED_SK_N_BE2c_H0FAN_L002',
        'DEDUP_SORTED_SK_N_FI_H0FAV_L007')

abbrev_ID<-c('CHLA90_L008',
               'COG16_L004',
               'COG21_L005',
               'COG25_L006',
               'COG28_L007',
               'COG29_L008',
               'COG7_L003',
               'LAN_6_MUT_P375_1',
               'LAN_6_pre_MUT_P81_2',
               'LAN_6_wt_P75_3',
               'LAN_6_wt_P259_4',
               'CHLA_291_C3_5',
               'CHLA_291_C4_6',
               'CHW_11_7',
               'COG_4_8',
               'COG_32_1',
               'LAN_6_MUT_P375_2',
               'LAN_6_pre_MUT_P81_3',
               'LAN_6_wt_P75_4',
               'CHLA_291_C3_2',
               'CHLA_291_C4_3',
               'CHW_11_4',
               'COG_4_5',
               'COG_32_6',
               'LAN_6_wt_P259_7',
               'SH_SY5Y_L001',
               'SK_N_BE2c_L002',
               'SK_N_FI_L007')

source_site<-c('Macrogen',
               'Macrogen',
               'Macrogen',
               'Macrogen',
               'Macrogen',
               'Macrogen',
               'Macrogen',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Garvan',
               'Macrogen',
               'Macrogen',
               'Macrogen')

db_source<-data.frame(source_key=keys, abbrev_ID, source_site=source_site)

# read trimmomatic log files

#Input Read Pairs: 360131678 Both Surviving: 272904363 (75.78%) Forward Only Surviving: 28908510 (8.03%) Reverse Only Surviving: 34119431 (9.47%) Dropped: 24199374 (6.72%)

input_dir<-'/mnt/TempSAN/eteber/TRIMMOMATIC_LOGS'
script<-c()
temp1<-list.files('/mnt/TempSAN/eteber/TRIMMOMATIC_LOGS', pattern=glob2rx('*.log'))
temp1<-gsub('.log','',temp1)

temp1<-unique(sort(temp1))
length(temp1)
sample_list<-temp1

library(stringr) 

db_trimmomatic<-data.frame()

for (sampleID in sample_list)
{
  
  sample_file     <-paste(sampleID,'.log',sep='')
  sh_string<-paste('cat', file.path(input_dir,sample_file), '| grep \"Input Read Pairs\"', sep=' ') 
  out<-system(sh_string, intern = TRUE)
  if (length(out)> 0)
  { 
      trimmed_line              <- str_trim(gsub("\\([^()]*\\)", "", out))
      array                     <- unlist(str_split(trimmed_line,' '))
      input_read_pairs          <-as.numeric(as.character(array[4]))
      both_surviving            <-as.numeric(as.character(array[7]))
      PCT_both_surviving        <-both_surviving/input_read_pairs 
      
      forward_only_surviving       <-as.numeric(as.character(array[12]))
      PCT_forward_only_surviving   <-forward_only_surviving/input_read_pairs 
      
      reverse_only_surviving        <-as.numeric(as.character(array[17]))
      PCT_reverse_only_surviving    <-reverse_only_surviving/input_read_pairs
      
      dropped                   <-as.numeric(as.character(array[20]))
      PCT_dropped               <-dropped/input_read_pairs 
    
      sample_key_trimmmomatic<-paste('DEDUP_SORT_', sampleID, sep='')
      trim_rec<-data.frame(sample_key_trimmmomatic, input_read_pairs,both_surviving, PCT_both_surviving, forward_only_surviving, PCT_forward_only_surviving, reverse_only_surviving, PCT_reverse_only_surviving,dropped,PCT_dropped)
      db_trimmomatic<-rbind(db_trimmomatic,trim_rec)  
  }
}


############################################################################


db_picard_wgstools
db_picard_alignment
db_picard_insert_size
db_samtools
db_fastqc


combined_table<-cbind(db_source,db_trimmomatic,db_picard_wgstools, db_picard_alignment, db_picard_insert_size, db_samtools, db_fastqc)

#combined_table<-sqldf('SELECT * FROM db_picard_wgstools INNER JOIN db_samtools ON TRIM(db_picard_wgstools.sample_name) = TRIM(db_samtools.sample_name)')

write.table(combined_table,file.path('/mnt/TempSAN/eteber/Sequencing_Mapping_Stats','read_mapping_stats.txt'), sep='\t')




pdfPath = '/mnt/san2/bioinformatics/eteber/cancer_group/rdagg/plots_figures/sequencing_mapping_quality.pdf'
#pdf(file=pdfPath, width=14, height=8) 

pdf(file=pdfPath, width=22, height=14) 

#par(mfrow=c(5,1)

## TRIMMOMATIC FIGURE

#sets the bottom, left, top and right margins respectively


#par(mar=c(10, 4 , 1, 12))
oma=c(1,1,1,1)
trim_report<-cbind(db_source, db_trimmomatic)
trim_report<-trim_report[order(source_site, abbrev_ID),]
stacker<-c('both_surviving','forward_only_surviving', 'reverse_only_surviving', 'dropped')
# Stacked Bar Plot with Colors and Legend
rownames(trim_report)<-trim_report$abbrev_ID
#par(mar=c(12,5,4,4),xpd=TRUE)
par(mar=c(12,5,4,14),xpd=TRUE)
counts<-t(trim_report[,stacker])
barplot(counts, main="after adapter removal and trimming reads", xlab="", col=1:length(stacker), las=2)
#y-axis = read_pairs
# Add legend to top right, outside plot region
legend("topright", inset=c(-0.1,0), legend=rownames(counts), title="read pairs", fill=c(1:length(stacker)))
#par(mar=c(1, 1 , 1, 1))
#plot.new()
#legend("topleft", legend=rownames(counts), title="read pairs", fill=c(1:length(stacker)))


#######################################################################

#filter<-c("PCT_EXC_MAPQ", "PCT_EXC_DUPE", "PCT_EXC_UNPAIRED", "PCT_EXC_BASEQ",  "PCT_EXC_OVERLAP",  "PCT_EXC_CAPPED", "PCT_EXC_TOTAL")
stacker<-c("PCT_EXC_MAPQ", "PCT_EXC_DUPE", "PCT_EXC_UNPAIRED", "PCT_EXC_BASEQ",  "PCT_EXC_OVERLAP",  "PCT_EXC_CAPPED")


wgs_report<-cbind(db_source,db_picard_wgstools[, stacker])
wgs_report<-wgs_report[order(source_site, abbrev_ID),]

rownames(wgs_report)<-wgs_report$abbrev_ID
par(mar=c(12,5,4,14),xpd=TRUE)
counts<-t(wgs_report[,stacker])
barplot(counts, main="fraction of aligned bases that were filtered out", xlab="", col=1:length(stacker), las=2)
#y-axis = read_pairs
# Add legend to top right, outside plot region
legend("topright", inset=c(-0.1,0), legend=rownames(counts), title="filtered out bases", fill=c(1:length(stacker)))


######################################################################

#coverage plot

var<-c("MEAN_COVERAGE", "SD_COVERAGE", "MEDIAN_COVERAGE", "MAD_COVERAGE")

wgs_report<-cbind(db_source,db_picard_wgstools[, var])
wgs_report<-wgs_report[order(source_site, abbrev_ID),]
rownames(wgs_report)<-wgs_report$abbrev_ID
x=1:length(rownames(wgs_report))
par(mar=c(12,5,4,14),xpd=TRUE)
sd.up = as.numeric(wgs_report$MEAN_COVERAGE)+as.numeric(wgs_report$SD_COVERAGE)
sd.dn = as.numeric(wgs_report$MEAN_COVERAGE)-as.numeric(wgs_report$SD_COVERAGE)
plot(wgs_report$MEAN_COVERAGE~x, cex=1.5, xaxt='n', xlab='',ylab='MEAN COVERAGE (STDEV)', col='black',pch=16, ylim=c(0,60))
axis(1, at=x, labels=rownames(wgs_report), las=2)
arrows(x,sd.dn,x,sd.up,code=3,length=0.1,angle=90,lty=1, lwd=2, col='red')


###############################################################################################

stacker <-c('total_reads',  'duplicate_reads',	'mapped_reads',	'all_paired_reads',	'count_reads1',	'count_reads2',	'proper_pair',	'with_itself_and_mate_mapped',	'singletons',	'with_mate_mapped_to_a_different_chr',	'with_mate_mapped_to_a_different_chr_mapQ5')

samtool_report<-cbind(db_source,db_samtools[, stacker])
samtool_report<-samtool_report[order(source_site, abbrev_ID),]

rownames(samtool_report)<-samtool_report$abbrev_ID
par(mar=c(12,5,4,14),xpd=TRUE)
counts<-t(samtool_report[,stacker])
barplot(counts, main="samtools", xlab="", col=1:length(stacker), las=2)
#y-axis = read_pairs
# Add legend to top right, outside plot region
legend("topright", inset=c(-0.15,0), legend=rownames(counts), title="samtools", cex=0.7, fill=c(1:length(stacker)))


#############################

var<-c('abbrev_ID','MEAN_READ_LENGTH', 'MEAN_INSERT_SIZE', 'STANDARD_DEVIATION', 'MEDIAN_INSERT_SIZE')

size_report<-combined_table[, var]
size_report<-size_report[order(source_site, abbrev_ID),]
rownames(size_report)<-size_report$abbrev_ID
x=1:length(rownames(size_report))
par(mar=c(12,5,4,14),xpd=TRUE)
sd.up = as.numeric(size_report$MEAN_INSERT_SIZE)+as.numeric(size_report$STANDARD_DEVIATION)
sd.dn = as.numeric(size_report$MEAN_INSERT_SIZE)-as.numeric(size_report$STANDARD_DEVIATION)
plot(size_report$MEAN_INSERT_SIZE~x, cex=1.5, xaxt='n', xlab='',ylab='READ LENGTH AND INSERT SIZE', col='black',pch=16, ylim=c(0,600))
points(size_report$MEAN_READ_LENGTH~x, cex=1.5, pch=21,  bg="lightgreen")
axis(1, at=x, labels=rownames(size_report), las=2)
arrows(x,sd.dn,x,sd.up,code=3,length=0.1,angle=90,lty=1, lwd=2, col='red')
legend("topright", inset=c(-0.15,0), legend=c('MEAN_INSERT_SIZE (STDEV)', 'MEAN_READ_LENGTH'), title="SIZE", cex=0.8, pt.cex=1.2,  pch=c(16,16), col=c('black', 'lightgreen'))



#plot.new()
#legend('bottomleft', legend=legend_numbers, pch = legend_numbers, col= legend_numbers+1, cex=0.85, horiz = TRUE, title = name)

dev.off() 
  
