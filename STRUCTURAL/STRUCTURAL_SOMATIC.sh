#!/bin/bash 

IX='/home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/8'

for file in $IX/SOMATIC_GRIDSS_*.vcf; do
    Rscript STRUCTURAL/STRUCTURAL_SOMATIC.R $file $file.bed
done
