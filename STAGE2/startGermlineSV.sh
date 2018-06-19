#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    x=$(qsub -v sampleID=$name STAGE2/germlineSV.pbs)
    echo $name $x
done;
