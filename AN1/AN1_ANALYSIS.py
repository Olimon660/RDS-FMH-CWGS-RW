#
# Python script for analyzinggermline resulting after filtering
#
#   python3 AN1/AN1_ANALYSIS.py AN1/FILTERED_AN1.csv > ANALYSIS_FILTERED_AN1.csv
#
# Generate a user-friendly TSV output file for the results.
#
#   - For each contrast result
#   - For each mortal and immortal sample
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

def analyze(cons, g, p):
    for c in cons:
        m1 = c["m1"] # WGS germline mortal
        m2 = c["m2"] # WGS germline immortal

        def check(m):
            # Only the germlines for the sample
            x = [i for i in g if i["samp"] == m]
            
            # No mutation? This is (+)
            if len(x) == 0:
                print(m + "+")
                
            # Mutation?
            else:
                print(m + "-")
        
        check(m1)
        check(m2)

#
# Parse WGS results and convert it into a format more useful for further analysis
#

def WGS(cons, line):
    toks = line.strip().split(';')        
        
    # Eg: 6/ANNOTATED_REMOVED_FILTERED_INDEL_NORM_DECOM_GATK_IIICF-T_B3.vcf
    samp = file2Samp(toks[1])

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

    return { "samp": samp, "chr":chr, "pos":pos, "ref":ref, "alt":alt, "imp":imp, "gene":gene }

#
# Parse proteomics data and convert it into a format more useful for further analysis
#

def PROTS(cons, line):
    pass

with open(sys.argv[1], "r") as r1:
    g = [] # WGS germline
    p = [] # Proteomics
    
    for l1 in r1:
        if len(l1.strip()) == 0:
            continue
        
        # Contrasts
        cons = []
        
        with open("AN1/AN1_CONTRASTS.csv") as r2:
            for l2 in r2:
                if "Mortal" in l2 and "Immortal" in l2:
                    continue
                toks = l2.strip().split(',')                
                cons.append({
                    "m1": toks[0],
                    "m2": toks[1]
                })
                
        # Search for WGS data
        g.append(WGS(cons, l1))

    analyze(cons, g, p)
