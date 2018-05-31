#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    #GATK=$(qsub -v sampleID=$name HAPLOTYPECALLER.pbs)
    echo $t_name $GATK
done;
