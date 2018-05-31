#!/bin/bash 

my_file='/project/RDS-SMS-WGS22-RW/haplotype_input_file_names.txt'
END=$(cat $my_file | wc -l)
START=1

for i in $(eval echo "{$START..$END}");
 do
   	t_name=$(head -"$i" $my_file | tail -1 | cut -f1)
   	echo $t_name
	
   	GATK=$(qsub -v sampleID=$t_name HAPLOTYPECALLER.pbs)
        echo $t_name $GATK
done;
