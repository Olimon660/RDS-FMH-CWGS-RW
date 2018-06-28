#!/bin/bash 

file=$1
samples=$(cat $file)

for name in 6/ANNOTATED_REMOVED_FILTERED_/*; do
    $(qsub -v file=$name GERMLINE/GERMLINE_DATABASE.pbs)
done;