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

def isMod(x):
    return x == "MODIFIER"
    
def isHigh(x):
    return x == "HIGH"

with open(sys.argv[1], "r") as r:
    for line in r:
        if len(line.strip()) == 0:
            continue
        toks = line.strip().split(';')        
        
        # Eg: 6/ANNOTATED_REMOVED_FILTERED_INDEL_NORM_DECOM_GATK_IIICF-T_B3.vcf
        file = toks[1]
        
        # Ensembl gene name
        gene = toks[19]
        
        # Gene symbol
        sym = toks[18]
        
        # Feature type
        ft = toks[20]
        
        # Feature
        ff = toks[21]
        
        # Impact (LOW, HIGH, MODERATE, MODIFIER)
        impact = toks[17]

        chr = toks[11]
        pos = toks[12]
        ref = toks[13]
        alt = toks[14]

        # Consequence
        con = toks[16]

        print(impact)
        

        asadsds
