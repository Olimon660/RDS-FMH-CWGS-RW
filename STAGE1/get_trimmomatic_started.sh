#
# ./STAGE1/get_trimmomatic_started.sh samples.txt
#

#!/bin/bash 

file=$1
samples=$(cat $file)

for name in $samples; do
    A=$(qsub -v sm=$name STAGE1/1_HPC_hg38_TRIMMOMATIC.pbs) 
    echo $name $A
done;
