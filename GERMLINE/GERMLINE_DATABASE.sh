#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -P RDS-FMH-CWGS-RW
#PBS -l select=1:ncpus=1:mem=32GB
#PBS -q defaultQ

module load python
cd /scratch/RDS-FMH-CWGS-RW/Sources

#python OTHERS/ANNOTATION_SQL.py OTHERS/ANNOTATION.sql /scratch/RDS-FMH-CWGS-RW/Sources/6 /scratch/RDS-FMH-CWGS-RW/Sources/6/germline.bin
