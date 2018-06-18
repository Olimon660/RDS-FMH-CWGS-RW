#!/bin/bash

NORMAL=REALIGNED_RG_DEDUP_SORTED_HG19_JFCF6_P12_H06L4ALXX_6.bam
TUMOUR=REALIGNED_RG_DEDUP_SORTED_HG19_JFCF6_T_1J_1-3C_H06L4ALXX_3.bam # The one Erdhal asked for

java -jar gridss-1.7.2-gridss-jar-with-dependencies.jar \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.compression_level=1 \WORKER_THREADS=15 \
	ASSEMBLY=gridss_JFCF6_SV.bam \
	INPUT="$NORMAL" \
	INPUT="$TUMOUR" \
	OUTPUT=gridss_JFCF6_SV.vcf \
	REFERENCE_SEQUENCE=combined.fasta &> gridss_JFCF6_SV.log
