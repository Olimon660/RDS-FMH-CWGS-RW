#
# ./STAGE1/2.sh <INPUT> <SAMPLE_TXT>
#

#!/bin/bash 

samples=$(cat $2)

for name in $samples; do
    print('qsub -v sampleID=$name -v in=$1 STAGE1/1_HPC_hg38_TRIMMOMATIC.pbs')
    qsub -v sampleID=$name -v in=$1 STAGE1/1_HPC_hg38_TRIMMOMATIC.pbs
done;
