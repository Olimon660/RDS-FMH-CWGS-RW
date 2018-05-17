#
# ./STAGE1/get_BWA.started.sh
#

#!/bin/bash 

samples=$(cat samples.txt)

for name in $samples; do
    A=$(qsub -v sm=$name STAGE1/4_HPC_hg38_alignment.pbs) 
    echo $name $A
done;
