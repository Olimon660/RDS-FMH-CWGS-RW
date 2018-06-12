#!/bin/bash 

my_file='/project/RDS-SMS-WGS22-RW/case_control_input_files.txt'
END=$(cat $my_file | wc -l)
START=1

for i in $(eval echo "{$START..$END}");
 do
   	t_name=$(head -"$i" $my_file | tail -1 | cut -f1)
   	n_name=$(head -"$i" $my_file | tail -1 | cut -f2)
   	echo $t_name $n_name

	if [ $n_name != "NONE" ]
	then
   		MUTECT=$(qsub -v sampleID=$t_name,controlID=$n_name MUTECT_with_control.pbs)
        	echo $t_name $n_name $MUTECT
	else
   		MUTECT=$(qsub -v sampleID=$t_name MUTECT_without_control.pbs)
        	echo $t_name $MUTECT
	fi
done;
