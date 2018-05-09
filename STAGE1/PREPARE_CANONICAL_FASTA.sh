hg38='/home/CMRI/eteber/das/Process/Bioinformatics/eteber/hg38bundle/hg38bundle/broad_38_canonical.fasta'
REF='/home/CMRI/eteber/das/Process/Bioinformatics/eteber/hg38bundle/hg38bundle/broad_hg38_canonical_with_vectors.fasta'
DS='/home/CMRI/eteber/das/Process/Bioinformatics/eteber/hg38bundle/hg38bundle'

SV40='pRSVT_seq.fasta'
HPVE6E7='HPV_E6E7_seq.fasta'
EBV='NC_007605_EBV.fasta'
pLKO1='pLKO1_scrambled_seq.fasta'
pLXSN='pLXSN_seq.fasta'
pBABE='pBABE_puro_seq.fasta'

cd $DS

#/tools/samtools/1.3/samtools faidx broad_38.fasta
#/tools/samtools/1.3/samtools faidx broad_38.fasta chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM  > broad_38_canonical.fasta

time cat $hg38 $EBV $SV40 $pLKO1 $pLXSN $HPVE6E7 $pBABE > $REF
 
SM='/tools/samtools/1.3'
MP='/tools/bwa/0.7.10'

time $SM/samtools faidx $REF 
time $MP/bwa index -a bwtsw $REF

cat $REF | grep '>'

###reference ready to use
