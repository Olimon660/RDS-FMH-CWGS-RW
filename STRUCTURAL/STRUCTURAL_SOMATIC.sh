#!/bin/bash

IS='/home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/5'
OX='/home/twong/das_manual/Process/Bioinformatics/tedwong/RDS-FMH-CWGS-RW/8'

NORM="$1"
TUMO="$2"

NORMAL=${IS}/REALIGNED_RG_DEDUP_SORTED_HG19_$NORM.bam
TUMOR=${IS}/REALIGNED_RG_DEDUP_SORTED_HG19_$TUMO.bam

echo $NORMAL
echo $TUMOR

OUTPUT=${OX}/SOMATIC_GRIDSS_$TUMO.vcf
ASSEMBLY=${OX}/SOMATIC_GRIDSS_$TUMO.bam

echo 'java -ea -Xmx128g \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.compression_level=1 \
	-cp gridss-1.7.2-gridss-jar-with-dependencies.jar gridss.CallVariants \
	TMP_DIR=8 \
	WORKER_THREADS=2 \
	WORKING_DIR=8 \
	REFERENCE_SEQUENCE=combined.fasta \
	INPUT=$NORMAL \
	INPUT=$TUMOR \
	OUTPUT=$OUTPUT \
	ASSEMBLY=$ASSEMBLY \
	2>&1 | tee -a gridss.somatic.$TUMO.$$.log'
