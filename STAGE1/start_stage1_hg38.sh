#!/bin/bash 

file=$1
samples=$(cat $file)

cd /project/RDS-FMH-CWGS-RW

for sname in $samples; do
	# MARK DUPS
	G=$(qsub -W depend=afterok:$F -v sampleID=$sname 7_HPC_hg38_markdups.pbs) 
	echo $sname $G 
	
	# DEL SORTED BAM
	H=$(qsub -W depend=afterok:$G -v sampleID=$sname 8_HPC_hg38_del_sorted_bam.pbs)
        echo $sname $H

	# READGROUPS
	I=$(qsub -W depend=afterok:$H -v sampleID=$sname 9_HPC_hg38_readgroups.pbs) 
	echo $sname $I

	# DEL_DEDUP
	J=$(qsub -W depend=afterok:$I -v sampleID=$sname 10_HPC_hg38_delete_dedup.pbs)
        echo $sname $J

	# INDEX
	K=$(qsub -W depend=afterok:$J -v sampleID=$sname 11_HPC_hg38_index.pbs) 
	echo $sname $K
 
	# RealignerTargetCreator
	L=$(qsub -W depend=afterok:$K -v sampleID=$sname 12_HPC_hg38_RealignerTarget.pbs) 
	echo $sname $L

	# IndelRealigner
	M=$(qsub -W depend=afterok:$L -v sampleID=$sname 13_HPC_hg38_IndelRealign.pbs) 
	echo $sname $M

        # FINAL INDEX
        N=$(qsub -W depend=afterok:$M -v sampleID=$sname 14_HPC_hg38_final_index.pbs)
        echo $sname $N
done;
