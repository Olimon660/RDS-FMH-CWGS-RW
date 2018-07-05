library(plyr)
library(GenomicFeatures)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#seqlevels(txdb)

##################################################################

#generate transcript list by gene - GRange
#gene_intervals <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by= "gene")

#gr_gene_intervals<-unlist(gene_intervals)
#gr_gene_intervals$GENEID<-names(gr_gene_intervals)

#length(gr_gene_intervals)

#############

#df_gene_transcripts<-data.frame(TXCHR=seqnames(gr_gene_intervals), ranges(gr_gene_intervals), strand=strand(gr_gene_intervals), values(gr_gene_intervals))
#colnames(df_gene_transcripts)[which(names(df_gene_transcripts) == "start")]   <- "TXSTART"
#colnames(df_gene_transcripts)[which(names(df_gene_transcripts) == "end")]     <- "TXEND"
#colnames(df_gene_transcripts)[which(names(df_gene_transcripts) == "strand")]  <- "TXSTRAND"
#colnames(df_gene_transcripts)[which(names(df_gene_transcripts) == "width")]   <-  "TXWIDTH"
#colnames(df_gene_transcripts)[which(names(df_gene_transcripts) == "GENEID")]  <-  "TXGENEID"
#colnames(df_gene_transcripts)[which(names(df_gene_transcripts) == "tx_name")] <-  "TX_NAME"
#colnames(df_gene_transcripts)[which(names(df_gene_transcripts) == "tx_id")]   <-  "TX_ID"

#clean up

#DROP <- c("names")

#KEEP <- !(colnames(df_gene_transcripts) %in% DROP)

#df_gene_transcripts<-df_gene_transcripts[,KEEP]

#head(df_gene_transcripts)


#library(org.Hs.eg.db)
#library(annotate)

#e2s = toTable(org.Hs.egSYMBOL)

# extract the longest transcript per gene
#library(plyr)

#db_longest_transcripts<-ddply(df_gene_transcripts,~TXGENEID,function(x){x[which.max(x$TXWIDTH),]})
#head(db_longest_transcripts)
#dim(db_longest_transcripts)

#clean up
#DROP <- c("names")
#KEEP <- !(colnames(db_longest_transcripts) %in% DROP)
#db_longest_transcripts<-db_longest_transcripts[,KEEP]

# create a Grange from data.frame
#longest_transcripts_txDB<-makeGRangesFromDataFrame(db_longest_transcripts, keep.extra.columns=TRUE)

#outDir<-'/processing/PROCESSING/MIGRATED/Bioinformatics/eteber/temp'
#dataDir<-'/processing/PROCESSING/MIGRATED/Bioinformatics/eteber/CN_24_SAMPLES_WGS/FN_DELLY_VCFS'


my_class_A_results <- function(x, txdb, df_gene_transcripts, e2s, outDir)
{
    vcf <- try(readVcf(x, "hg19"))
    vcf  <-vcf[fixed(vcf)$FILTER =='PASS',]

    info_columns<-rownames(info(header(vcf)))

    vchr1                        <-as.character(seqnames(vcf))
  vchr2                        <-info(vcf)['CHR2']$CHR2
  vstart                       <-start(vcf)
  vend                         <-info(vcf)['END']$END
  vSVTYPE                      <-info(vcf)['SVTYPE']$SVTYPE

  ##############################################

  vID                          <-names(rowRanges(vcf))
  vQUAL                        <-qual(vcf)
  vREF                         <-as.character(ref(vcf))
  vALT                         <-as.character(alt(vcf))
  #ALT                         <-CharacterList(ALT)
  vconnection_type             <-info(vcf)['CT']$CT
  vread_pairs_support          <-info(vcf)['PE']$PE
  vsoft_reads_support          <- info(vcf)['SR']$SR
  vsoft_read_consensus_aln_qty <- info(vcf)['SRQ']$SRQ
  vmap_quality                 <- info(vcf)['MAPQ']$MAPQ
  pre                         <- info(vcf)['PRECISE']
  im                          <- info(vcf)['IMPRECISE']
  vresolution                  <- ifelse(pre$PRECISE==TRUE, 'PRECISE', 'IMPRECISE')
  diff                        <-vend - vstart + 1
  vsize                        <- ifelse(vSVTYPE=='TRA',NA,diff)

  # when there is a control DELLY prodcues an RDRATIO

  my_function <-function(x, colnames) {
    rdratio_idx<-which(colnames=='RDRATIO')
    if (!is.null(rdratio_idx)) {return(info(x)[['RDRATIO']])}
    else {return}
  }

  rdratio_idx<-which(info_columns=='RDRATIO')

  if (length(rdratio_idx) > 0) {
    vRDRATIO<-info(vcf)[['RDRATIO']]
  } else  { vRDRATIO<-rep(NA, length(vcf))}

  db_vcf_entries<-data.frame(chr1=vchr1, start=vstart, chr2=vchr2, end=vend, SVID=vID, SVSIZE = vsize, REF = vREF, ALT = vALT, SVTYPE = vSVTYPE, resolution  = vresolution , connection_type = vconnection_type, read_pairs_support =vread_pairs_support, soft_reads_support =vsoft_reads_support, soft_read_consensus_aln_qty=vsoft_read_consensus_aln_qty,  median_map_quality = vmap_quality, RDRATIO=vRDRATIO)

  first_breakpoint <-with(db_vcf_entries, GRanges(seqnames = chr1,
                                                  ranges = IRanges(start=start, end=start),
                                                  strand = NULL,
                                                  SVID,
                                                  SVSIZE,
                                                  REF,
                                                  ALT,
                                                  SVTYPE,
                                                  resolution,
                                                  connection_type,
                                                  read_pairs_support,
                                                  soft_reads_support,
                                                  soft_read_consensus_aln_qty,
                                                  median_map_quality,
                                                  RDRATIO))

  second_breakpoint <- with(db_vcf_entries, GRanges(seqnames = chr2,
                                                    ranges = IRanges(start=end, end=end),
                                                    strand = NULL,
                                                    SVID,
                                                    SVSIZE,
                                                    REF,
                                                    ALT,
                                                    SVTYPE,
                                                    resolution,
                                                    connection_type,
                                                    read_pairs_support,
                                                    soft_reads_support,
                                                    soft_read_consensus_aln_qty,
                                                    median_map_quality,
                                                    RDRATIO))

  #set up CLASS_A GRANGES OBJECT
  tmp<-GRangesList(first_breakpoint, second_breakpoint)

  CLASS_A<-unlist(tmp)

  db_CLASS_A<-data.frame(SVCHR=seqnames(CLASS_A), ranges(CLASS_A), mcols(CLASS_A))
  db_CLASS_A$QUERY_NO<-c(1:length(CLASS_A))

  colnames(db_CLASS_A)[which(names(db_CLASS_A) == "start")] <- "SVSTART"
  colnames(db_CLASS_A)[which(names(db_CLASS_A) == "end")]   <- "SVEND"

  DROP <- c("width")
  KEEP <- !(colnames(db_CLASS_A) %in% DROP)

  db_CLASS_A<-db_CLASS_A[,KEEP]
  AllVariants(intergenic=IntergenicVariants(2000, 2000))

  

  
  # locate all variants using CLAS_A breakpoints
  HITS_A <- locateVariants(CLASS_A, txdb, AllVariants())

  # add SV_ID column to HITS GRange

  SVID<-values(CLASS_A[HITS_A$QUERYID, 'SVID'])$SVID

  HITS_A$HITSVID<-SVID

 

  # makie the HITS into a data.frame

 

  db_HITS_A<-data.frame(ranges(HITS_A), strand=strand(HITS_A), mcols(HITS_A))

  colnames(db_HITS_A)[which(names(db_HITS_A) == "start")]  <- "HITSTART"

  colnames(db_HITS_A)[which(names(db_HITS_A) == "end")]    <- "HITEND"

  colnames(db_HITS_A)[which(names(db_HITS_A) == "strand")] <- "HITSTRAND"

  colnames(db_HITS_A)[which(names(db_HITS_A) == "width")]  <- "HITWIDTH"

  colnames(db_HITS_A)[which(names(db_HITS_A) == "GENEID")]  <- "HITGENEID"

 

  head(db_HITS_A)

  head(db_CLASS_A)

 

  ######

 

  #bln_gne_of_interest<-db_HITS_A$HITGENEID==144455

  #db_HITS_A[which(bln_gne_of_interest==TRUE),]

 

  

  # MERGE HITS WITH CLASS_A DETAILS

 

  temp<-merge(db_HITS_A, db_CLASS_A, by.x='QUERYID', by.y='QUERY_NO')

  head(temp)

 

  #head(df_gene_transcripts)

  # MERGE HITS WITH GENE DETIALS DETAILS

 

  db_SV_RESULTS.CLASS_A<-merge(temp, df_gene_transcripts, by.x='TXID', by.y='TX_ID',all.x=TRUE, sort=FALSE)

  head(db_SV_RESULTS.CLASS_A)

  dim(db_SV_RESULTS.CLASS_A)

 

  # OBTAIN THE SYMBOL NAME

 

   entrezGenes <- db_SV_RESULTS.CLASS_A$TXGENEID

  tb_entrezID2genename<-select(org.Hs.eg.db, keys = entrezGenes, keytype = "ENTREZID", columns = c("SYMBOL","ENTREZID"))

  db_SV_RESULTS.CLASS_A$SYMBOL<-tb_entrezID2genename$SYMBOL

 

  # CLEAN UP DATA>FRAME

 

  DROP<-c('QUERYID',  'HITSTART',    'HITEND', 'HITWIDTH', 'HITSTRAND', 'HITSVID', 'HITGENEID')

  KEEP <- !(colnames(db_SV_RESULTS.CLASS_A) %in% DROP)

  db_SV_RESULTS.CLASS_A<-db_SV_RESULTS.CLASS_A[,KEEP]

  head(db_SV_RESULTS.CLASS_A)

 

 


  PRECEDE<-ldply(db_SV_RESULTS.CLASS_A$PRECEDEID, toString)

  db_SV_RESULTS.CLASS_A$PRECEDE_GENEID<-PRECEDE$V1

 

  FOLLOW<-ldply(db_SV_RESULTS.CLASS_A$FOLLOWID, toString)

  db_SV_RESULTS.CLASS_A$FOLLOW_GENEID<-FOLLOW$V1

 

  select_cols <- c("SVTYPE", "SVID", "SVCHR", "SVSTART", "SVEND", "SVSIZE", "SYMBOL", "TXGENEID", "LOCATION", "TXID", "TX_NAME","TXCHR", "TXSTART", "TXEND", "TXWIDTH", "TXSTRAND",  "resolution", "connection_type", "read_pairs_support", "soft_reads_support", "soft_read_consensus_aln_qty", "median_map_quality", "RDRATIO", "PRECEDE_GENEID","FOLLOW_GENEID")

 

  db_SV_RESULTS.CLASS_A<-db_SV_RESULTS.CLASS_A[ ,select_cols]

  

  outname<-gsub('.vcf', '',basename(x))

  output_filename<-file.path(outDir,paste('CLASS_A_ANNOTATIONS_', outname, '.txt', sep=''))
  write.table(db_SV_RESULTS.CLASS_A, output_filename , sep='\t', col.names = TRUE, row.names = FALSE, quote=FALSE)
}

dellyFileList <- list.files('/Users/twong/Sources/RDS-FMH-CWGS-RW/ERDAHL', full.names = T)

for (item in dellyFileList)
{
  print (item)
  print (gsub('.FILTERED.vcf', '',basename(item)))
  my_class_A_results(item, txdb, df_gene_transcripts, e2s, outDir)
}


head(df_gene_transcripts)
