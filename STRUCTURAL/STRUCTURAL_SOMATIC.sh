#!/bin/bash 

IX='/scratch/RDS-FMH-CWGS-RW/Sources/8'
OX='/scratch/RDS-FMH-CWGS-RW/Sources/8'

module load R

for file in $IX/SOMATIC_GRIDSS_*.vcf; do
    Rscript STRUCTURAL/STRUCTURAL_SOMATIC.R $file
done
