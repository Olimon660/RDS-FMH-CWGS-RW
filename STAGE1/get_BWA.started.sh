#
# ./STAGE1/get_BWA.started.sh samples1.txt
# ./STAGE1/get_BWA.started.sh samples2.txt
#

#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    A=$(qsub -v sm=$name STAGE1/4_HPC_hg38_alignment.pbs) 
    echo $name $A
done;
