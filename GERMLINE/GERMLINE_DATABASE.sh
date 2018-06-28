#!/bin/bash 

for file in 6/ANNOTATED_REMOVED_FILTERED_*vcf; do
    echo $file
#    $(qsub -v file=$file OTHERS/DATABASE.pbs)
done;
