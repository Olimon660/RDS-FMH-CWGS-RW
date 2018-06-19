#!/bin/bash 

OUT='/home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/6'

for file in $OUT/*.vcf; do
    echo vt decompose -s $file -o $OUT/DECOMPOSED_$(basename "$file")
    vt decompose -s $file -o $OUT/DECOMPOSED_$(basename "$file")
done
