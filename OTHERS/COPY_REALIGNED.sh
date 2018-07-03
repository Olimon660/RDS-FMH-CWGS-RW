#!/bin/bash 

OUT='/scratch/RDS-FMH-CWGS-RW/Sources/5'

for file in $OUT/REALIGNED_RG_DEDUP_SORTED_HG19_*.bam; do    
    #qsub -v file=$file OTHERS/COPY_REALIGNED.pbs
    echo 'qsub -v file="$file" OTHERS/COPY_REALIGNED.pbs'
done
