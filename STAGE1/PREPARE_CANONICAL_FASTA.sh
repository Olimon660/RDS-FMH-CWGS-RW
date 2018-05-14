#!/bin/bash

#PBS -l walltime=10:00:00
#PBS -P RDS-FMH-CWGS-RW
#PBS -l select=1:ncpus=1:mem=128GB
#PBS -q defaultQ

module load samtools

# Reference FASTA
FASTA='/scratch/RDS-FMH-CWGS-RW/Sources/combined.fasta'

# Where for BWA
BWA='/scratch/RDS-FMH-CWGS-RW/softwares/bwa-0.7.17/bwa'

# Where for SAMTools
SAMTOOLS='samtools'

#/tools/samtools/1.3/samtools faidx broad_38.fasta
#/tools/samtools/1.3/samtools faidx broad_38.fasta chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM  > broad_38_canonical.fasta

$SAMTOOLS faidx $FASTA 
$BWA index -a bwtsw $FASTA
