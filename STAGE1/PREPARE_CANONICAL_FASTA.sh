#!/bin/bash

#PBS -l walltime=10:00:00
#PBS -P RDS-FMH-CWGS-RW
#PBS -l select=1:ncpus=1
#PBS -q defaultQ

module load samtools

# Reference FASTA
FASTA='/home/CMRI/eteber/das/Process/Bioinformatics/eteber/hg38bundle/hg38bundle/broad_38_canonical.fasta'

# Where for BWA
BWA='/tools/bwa/0.7.10'

# Where for SAMTools
SAMTOOLS='samtools'

#/tools/samtools/1.3/samtools faidx broad_38.fasta
#/tools/samtools/1.3/samtools faidx broad_38.fasta chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM  > broad_38_canonical.fasta

time $SAMTOOLS faidx $FASTA 
time $BWA index -a bwtsw $FASTA
