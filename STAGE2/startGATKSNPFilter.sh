#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    GATK=$(qsub -v sampleID=$name STAGE2/HAPLOTYPECALLER.pbs)
    echo $name $GATK
done;
