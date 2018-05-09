## Summary

Sample usage for the software:

`java -jar trimmomatic-0.33.jar PE -threads 8 -phred33 $IX/${short_name}_R1.fq.gz $IX/${short_name}_R2.fq.gz \
$OX/F1_PAIRED_${short_name}.fq $OX/F1_UNPAIRED_${short_name}.fq \
$OX/F2_PAIRED_${short_name}.fq $OX/F2_UNPAIRED_${short_name}.fq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75 2> $OX/TRIMMOMATIC_${short_name}.log`
