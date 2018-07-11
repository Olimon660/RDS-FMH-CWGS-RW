#!/bin/bash 

IX='/scratch/RDS-FMH-CWGS-RW/Sources/5'

for file in $IX/REALIGNED*.bam; do
for file in $IX/REALIGNED_RG_DEDUP_SORTED_HG19_*ZK-58.bam; do
    qsub -v file=$file STATS/STATS_ALIGNMENT.pbs
done
