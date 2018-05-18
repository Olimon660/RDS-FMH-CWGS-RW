#
# ./STAGE1/get_BWA.started.sh
#

#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    A=$(qsub -v sm=$name STAGE1/4_HPC_hg38_alignment.pbs) 
    echo $name $A
done;
