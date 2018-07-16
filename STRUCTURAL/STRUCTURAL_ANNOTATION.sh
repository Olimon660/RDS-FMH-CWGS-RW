#!/bin/bash

IX='8'
OX='8'

for file in $IX/SOMATIC_GRIDSS_*.vcf; do
    java -Xmx16g -jar /home/twong/snpEff/snpEff.jar hg19 $IX/$(basename ${file}) > $OX/ANNOTATED_$(basename ${file})
done
