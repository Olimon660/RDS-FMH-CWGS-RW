#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    GATK=$(qsub -v sampleID=$name STAGE2/GATK_Filter_SNP.pbs)
    echo $name $GATK
done;
