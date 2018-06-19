#!/bin/bash 

file=$1
samples=$(cat $file)

for file in /home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/6/*.vcf; do
    echo $file
done
