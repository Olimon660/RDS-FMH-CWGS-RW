#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -P RDS-FMH-CWGS-RW
#PBS -l select=1:ncpus=1:mem=32GB
#PBS -q defaultQ

cd /scratch/RDS-FMH-CWGS-RW/Sources
module load python
python GERMLINE/ANNOTATION_SQL.py GERMLINE/ANNOTATION.sql 6 6/6.csv
