#!/bin/bash

#PBS -l walltime=12:00:00
#PBS -P RDS-FMH-CWGS-RW
#PBS -l select=1:ncpus=1:mem=32GB
#PBS -q small-express

module load python

cd /scratch/RDS-FMH-CWGS-RW/Sources
python OTHERS/ANNOTATION_SQL.py OTHERS/ANNOTATION.sql 7 7/7.csv
sqlite3 /scratch/RDS-FMH-CWGS-RW/Sources/7/7.db < /scratch/RDS-FMH-CWGS-RW/Sources/Somatic/SOMATIC_SQLLITE.sh
