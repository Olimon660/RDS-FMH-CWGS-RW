#
# Python script for analyzinggermline resulting after filtering
#
#   python3 AN1_ANALYSIS.py FILTERED_AN1.csv > ANALYSIS_FILTERED_AN1.csv
#

import sys

# Eg: 6/ANNOTATED_REMOVED_FILTERED_INDEL_NORM_DECOM_GATK_IIICF-T_B3.vcf
def file2Samp(x):
    assert("GATK" in x)
    return x.split("GATK_")[1].replace(".vcf", "")

with open(sys.argv[1], "r") as r:
    for line in r:
        toks = line.strip().split(';')
        
        print(line)
        
        #print(file2Samp(toks[1]))
        
        
