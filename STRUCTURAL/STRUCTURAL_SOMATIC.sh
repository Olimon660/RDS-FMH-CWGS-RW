#!/bin/bash 

IX='8'

for file in $IX/SOMATIC_GRIDSS_*.vcf; do
    Rscript STRUCTURAL/STRUCTURAL_SOMATIC.R $file $file.bed
done
