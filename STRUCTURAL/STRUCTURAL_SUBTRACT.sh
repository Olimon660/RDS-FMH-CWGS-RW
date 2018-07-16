#!/bin/bash 

IX='8'
OX='8'

for file in $IX/ANNOTATED_SOMATIC_GRIDSS_*.vcf; do
    Rscript STRUCTURAL/STRUCTURAL_SUBTRACT.R $file $OX/SUBTRACTED_$(basename ${file}).bed
done
