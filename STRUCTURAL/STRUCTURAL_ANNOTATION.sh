#!/bin/bash 

IX='/scratch/RDS-FMH-CWGS-RW/Sources/6'
OX='/scratch/RDS-FMH-CWGS-RW/Sources/6'

for file in $IX/REMOVED_FILTERED_*.vcf; do
    qsub -v file=$file GERMLINE/GERMLINE_ANNOTATION.pbs
done
