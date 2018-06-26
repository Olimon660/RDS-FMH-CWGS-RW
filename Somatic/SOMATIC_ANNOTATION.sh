#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    r=$(qsub -v file=$name Somatic/SOMATIC_ANNOTATION.pbs)
    echo $name $r
done;
