#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    t_name=$(head -"$i" $my_file | tail -1 | cut -f1)
    echo $t_name
	
    GATK=$(qsub -v sampleID=$t_name HAPLOTYPECALLER.pbs)
    echo $t_name $GATK
done;
