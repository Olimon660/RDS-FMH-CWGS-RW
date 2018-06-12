#!/bin/bash 

my_file='/scratch/RDS-FMH-CWGS-RW/Sources/STAGE2/case_control_input_files.txt'
END=$(cat $my_file | wc -l)
START=1

for i in $(eval echo "{$START..$END}");
    do
        normal=$(head -"$i" $my_file | tail -1 | cut -f1)
	tumor=$(head -"$i" $my_file | tail -1 | cut -f2)
	echo $normal $tumor

        MUTECT=$(qsub -v normalID=$normal,tumorID=$tumor MUTECT_with_control.pbs)
        echo $t_name $n_name $MUTECT
done;
