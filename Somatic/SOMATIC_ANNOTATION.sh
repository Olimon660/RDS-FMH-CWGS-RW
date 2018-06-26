#!/bin/bash 

IX='/scratch/RDS-FMH-CWGS-RW/Sources/7'
OX='/scratch/RDS-FMH-CWGS-RW/Sources/7'

for file in $IX/REMOVED_FILTERED_*.vcf; do
    qsub -v file=$file Somatic/SOMATIC_ANNOTATION.pbs
done
