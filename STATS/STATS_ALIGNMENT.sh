#!/bin/bash 

IX='/scratch/RDS-FMH-CWGS-RW/Sources/5'

for file in $IX/REALIGNED*.bam; do
    qsub -v file=$file STATS/STATS_ALIGNMENT.pbs
done
