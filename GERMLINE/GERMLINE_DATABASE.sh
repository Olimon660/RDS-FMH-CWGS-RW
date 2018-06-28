#!/bin/bash 

#for file in 6/ANNOTATED_REMOVED_FILTERED_*vcf; do
for file in 6/REMOVED_FILTERED_INDEL_NORM_DECOM_GATK_ZK-58*.vcf; do
    qsub -v file=$file OTHERS/DATABASE.pbs
done;
