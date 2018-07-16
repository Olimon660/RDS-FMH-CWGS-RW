#!/bin/bash

IX='/scratch/RDS-FMH-CWGS-RW/Sources/8'
OX='/scratch/RDS-FMH-CWGS-RW/Sources/8'

for file in $IX/SOMATIC_GRIDSS_*.vcf; do
    java -Xmx16g -jar /home/twong/snpEff/snpEff.jar hg19 $OX/$(basename ${file})
done
