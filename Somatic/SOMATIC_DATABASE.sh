#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -P RDS-FMH-CWGS-RW
#PBS -l select=1:ncpus=1:mem=32GB
#PBS -q defaultQ

module load python
python OTHERS/ANNOTATION_SQL.py OTHERS/ANNOTATION.sql 7 7/7.csv
