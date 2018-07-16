#!/bin/bash 

IX='/scratch/RDS-FMH-CWGS-RW/Sources/8'
OX='/scratch/RDS-FMH-CWGS-RW/Sources/8'

module load R

for file in $IX/SOMATIC_GRIDSS_*.vcf; do
    qsub -v file=$file STRUCTURAL/STRUCTURAL_SOMATIC.pbs
done
