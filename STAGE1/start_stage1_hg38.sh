#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    # MARK DUPS
    G=$(qsub -v sampleID=$name STAGE1/7_HPC_hg38_markdups.pbs) 
    echo $sname $G 
	
    # DEL SORTED BAM
    #H=$(qsub -W depend=afterok:$G -v sampleID=$name STAGE1/8_HPC_hg38_del_sorted_bam.pbs)
    #echo $sname $H

    # READGROUPS
    I=$(qsub -W depend=afterok:$G -v sampleID=$name STAGE1/9_HPC_hg38_readgroups.pbs) 
    echo $sname $I

    # DEL_DEDUP
    J=$(qsub -W depend=afterok:$I -v sampleID=$name STAGE1/10_HPC_hg38_delete_dedup.pbs)
    echo $sname $J
    
    # INDEX
    K=$(qsub -W depend=afterok:$J -v sampleID=$name STAGE1/11_HPC_hg38_index.pbs) 
    echo $sname $K
	
    # RealignerTargetCreator
    L=$(qsub -W depend=afterok:$K -v sampleID=$name STAGE1/12_HPC_hg38_RealignerTarget.pbs) 
    echo $sname $L

    # IndelRealigner
    M=$(qsub -W depend=afterok:$L -v sampleID=$name STAGE1/13_HPC_hg38_IndelRealign.pbs) 
    echo $sname $M

    # FINAL INDEX
    N=$(qsub -W depend=afterok:$M -v sampleID=$name STAGE1/14_HPC_hg38_final_index.pbs)
    echo $sname $N
done;
