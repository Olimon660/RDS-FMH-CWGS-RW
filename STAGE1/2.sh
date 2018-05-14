#
# ./STAGE1/2.sh <SAMPLE_TXT>
#

#!/bin/bash 

samples=$(cat $1)

for name in $samples; do
    A=$(qsub -v sm=$name STAGE1/1_HPC_hg38_TRIMMOMATIC.pbs) 
    echo $name $A
done;
