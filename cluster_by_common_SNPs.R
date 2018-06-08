library(VariantAnnotation)
library(stringr)


dir='/home/CMRI/eteber/WORKING_DIR/TEMP'

original<-list.files(dir, pattern=glob2rx('PASS_*.vcf'))
inputID<-gsub('.vcf','',original)
temp1<-gsub('PASS_SNV_GATK_REALIGNED_RG_DEDUP_SORTED_','',inputID)
sample_list<-temp1

df_ID<-data.frame(inputID=inputID,sampleID=sample_list)

ptm <- proc.time()

HASH_TABLE     <-data.frame()

for (sampleID in sample_list)
{
  
  filename<- paste(df_ID$inputID[df_ID$sampleID==sampleID], '.vcf', sep='')
  
  vcfFile = file.path(dir, filename) 
  vcf = readVcf(vcfFile, "hg19")
  
  CHR<-as.character(seqnames(vcf))
  POS<-start(ranges(vcf))
  REF<-as.character(ref(vcf))
  ALT<-as.character(unlist(alt(vcf)))
  
  ID<-paste(CHR,POS,REF,ALT,sep='_')
  
  SAMPLE_HASH<-data.frame(ID, sampleID)
  
  HASH_TABLE    <-rbind(HASH_TABLE,SAMPLE_HASH)
  print(sampleID)
}

proc.time() - ptm


#common SNPs

dbSNPfile<-'/home/CMRI/eteber/das/SAN/bioinformatics/eteber/hg19/dbsnp_138.left_normalised.hg19.sorted.COMMON.vcf.gz'
commonSNPfile = dbSNPfile
SNPvcf = readVcf(commonSNPfile, "hg19")

CHR<-as.character(seqnames(SNPvcf))
POS<-start(ranges(SNPvcf))
REF<-as.character(ref(SNPvcf))
ALT<-as.character(unlist(alt(SNPvcf)))

dbSNPID<-paste(CHR,POS,REF,ALT,sep='_')

length(dbSNPID)

match<-which(HASH_TABLE$ID %in% dbSNPID)
HASH_SNPs<-HASH_TABLE$ID[match]

ID<-sort(unique(HASH_SNPs))
length(ID)
no_unique_ID<-length(ID)
no_samples  <-length(sample_list)

# generate matrix
matrix_ID<-matrix(data=0, nrow=no_unique_ID, ncol=no_samples)
rownames(matrix_ID)<-ID
colnames(matrix_ID)<-sample_list

DB_SNP_SET<-HASH_TABLE[match,]

ix<-unique(cbind(as.character(DB_SNP_SET$ID),as.character(DB_SNP_SET$sampleID)))
matrix_ID[ix]<-1

#plot hierarchial cluster.
d = dist(t(matrix_ID), method = "binary")
hc = hclust(d, method="ward")
plot(hc,main='cluster based on common SNPs')

