#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    GATK=$(qsub -v sampleID=$name STAGE2/GERMLINE_FILTER_INDEL.pbs)
    echo $name $GATK
done;
