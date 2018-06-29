#!/bin/bash

DS='/home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/5'
OX='/home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/8'
REF='/home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/combined.fasta'

NORMAL=${DS}/REALIGNED_RG_DEDUP_SORTED_HG19_${normalID}.bam
TUMOR=${DS}/REALIGNED_RG_DEDUP_SORTED_HG19_${tumorID}.bam

REFERENCE=$REF
OUTPUT=${OX}/SOMATIC_GRIDSS_${tumorID}.vcf
ASSEMBLY=${OX}/SOMATIC_GRIDSS_${tumorID}.bam

echo $tumorID
cd /scratch/RDS-FMH-CWGS-RW/Sources/8

java -ea -Xmx128g \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.compression_level=1 \
	-cp /scratch/RDS-FMH-CWGS-RW/softwares/gridss-1.7.2-gridss-jar-with-dependencies.jar gridss.CallVariants \
	TMP_DIR=/scratch/RDS-FMH-CWGS-RW/Sources/8 \
	WORKER_THREADS=2 \
	WORKING_DIR=/scratch/RDS-FMH-CWGS-RW/Sources/8 \
	REFERENCE_SEQUENCE="$REFERENCE" \
	INPUT="$NORMAL" \
	INPUT="$TUMOR" \
	OUTPUT="$OUTPUT" \
	ASSEMBLY="$ASSEMBLY" \
	2>&1 | tee -a gridss.somatic.$HOSTNAME.$$.log
