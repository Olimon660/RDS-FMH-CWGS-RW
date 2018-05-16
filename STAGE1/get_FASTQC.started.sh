#
# ./STAGE1/get_FASTQC.started.sh
#

#!/bin/bash 

samples=$(cat samples.txt)

for name in $samples; do
    A=$(qsub -v sm=$name STAGE1/3_HPC_hg38_FASTQC.pbs) 
    echo $name $A
done;
