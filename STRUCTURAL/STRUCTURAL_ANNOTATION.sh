#!/bin/bash

IX='/home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/8'
OX='/home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/8'

for file in $IX/SOMATIC_GRIDSS_*.vcf; do
    java -Xmx16g -jar /home/twong/snpEff/snpEff.jar hg19 $IX/$(basename ${file}) > $OX/$(basename ${file}).txt
done
