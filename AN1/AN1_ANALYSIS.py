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
    return x == "MODERATE"
    
def isHigh(x):
    return x == "HIGH"
    
def isLow(x):
    return x == "LOW" or x == "MODIFIER"
    
def ens2Name(x):
    if x == "ENSG00000204209":
        return "DAXX"
    elif x == "ENSG00000085224":
        return "ATRX"
    elif x == "ENSG00000164362":
        return "TERT"
    elif x == "ENSG00000141510":
        return "TP53"

def WGS(line):
        toks = line.strip().split(';')        
        
        # Eg: 6/ANNOTATED_REMOVED_FILTERED_INDEL_NORM_DECOM_GATK_IIICF-T_B3.vcf
        file = toks[1]
        
        # Translated gene name (from Ensembl)
        gene = ens2Name(toks[19])

        # Gene symbol
        sym = toks[18]
        
        # Feature type
        ft = toks[20]
        
        # Feature
        ff = toks[21]
        
        # Impact (LOW, HIGH, MODERATE, MODIFIER)
        imp = toks[17]

        assert(isMod(imp) or isHigh(imp) or isLow(imp))

        chr = toks[11] # Chromosome
        pos = toks[12] # Position
        ref = toks[13] # Reference
        alt = toks[14] # Allele

        # Consequence
        con = toks[16]

with open(sys.argv[1], "r") as r:
    for line in r:
        if len(line.strip()) == 0:
            continue
        WGS(line)
