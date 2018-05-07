#!/bin/bash 
#cd /scratch/RDS-SMS-WGS22-RW/
directory=$1
#cd $directory
#filename=$(ls *_R1.fq.gz)
#remove substring
#my_sample="${filename//_R1.fq.gz}"
my_samples=$(cat $directory/sample_list_batch2.txt)
cd /project/RDS-SMS-WGS22-RW
for sname in $my_samples; do
	#TRIMMOMATIC
	A=$(qsub -v sampleID=$sname 1_HPC_hg38_TRIMMOMATIC.pbs) 
	echo $sname $A '1_HPC_hg38_TRIMMOMATIC'
done;
