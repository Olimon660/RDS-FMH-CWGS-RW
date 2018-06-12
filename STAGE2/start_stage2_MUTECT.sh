#!/bin/bash 

my_file='/scratch/RDS-FMH-CWGS-RW/Sources/STAGE2/case_control_input_files.txt'
END=$(cat $my_file | wc -l)
START=1

for i in $(eval echo "{$START..$END}");
    do
        t_name=$(head -"$i" $my_file | tail -1 | cut -f1)
	n_name=$(head -"$i" $my_file | tail -1 | cut -f2)
	echo $t_name $n_name

        MUTECT=$(qsub -v sampleID=$t_name,controlID=$n_name MUTECT_with_control.pbs)
        echo $t_name $n_name $MUTECT
done;
