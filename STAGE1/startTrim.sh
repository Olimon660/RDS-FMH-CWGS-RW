#
# ./STAGE1/startTrim.sh STAGE1
#

#!/bin/bash 

directory=$1
my_samples=$(cat $directory/samples.txt)

for sname in $my_samples; do
    A=$(qsub -v sampleID=$sname STAGE1/1_HPC_hg38_TRIMMOMATIC.pbs) 
done;
