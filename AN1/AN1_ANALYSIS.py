#
# Python script for analyzing resulting after filtering
#
#     python3 AN1/AN1_ANALYSIS.py F
#     python3 AN1/AN1_ANALYSIS.py W
#     python3 AN1/AN1_ANALYSIS.py P
#     python3 AN1/AN1_ANALYSIS.py A
#
# Generate a user-friendly TSV output file for the results.
#

import sys
import pickle

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

def load(x):
    with open(x, "rb") as r:
        return pickle.load(r)

def save(file, x):
    with open(file, "wb") as w:
        pickle.dump(x, w, protocol=pickle.HIGHEST_PROTOCOL)    

def WGSByName(x, name):
    return [i for i in x if i["name"] == name]

def onlySNP(x):
    return [i for i in x if i["type"] == "SNP"]

def onlyInd(x):
    return [i for i in x if i["type"] == "Ind"]

def onlyGene(x, gn):
    return [i for i in x if i["gn"] == gn]

def only(x, key, v):
    return [i for i in x if i[key] == v]

def analyze(W, P):
    # Format string
    f = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"

    # Genes interested
    genes = [ "ATRX", "DAXX" ]

    w = open("AN1/AN1_RESULTS.tsv", "w")
    w.write(f.format("Name", "Mortal", "Immortal", "Gene/Protein", \
                     "W_M_SNP", "W_M_Ind", "W_I_SNP", "W_I_Ind", "W_M_S", "W_I_S", \
                     "P_M_A", "P_I_A"))

    with open("AN1/AN1_CONTRASTS.csv") as r:
        for l in r:
            if "Mortal" in l and "Immortal" in l:
                continue
            toks = l.strip().split(',')

            m1  = toks[0] # Mortal
            m2  = toks[1] # Immortal
            m1W = WGSByName(W, m1)
            m2W = WGSByName(W, m2)
                
            # Name of the contrast
            name = m1 + "_" + m2

            W_M_SNP = onlySNP(m1W) # SNPs for mortal (all genes)
            W_M_IND = onlyInd(m1W) # Indels for mortal (all genes)
            W_I_SNP = onlySNP(m2W) # SNPs for immortal (all genes)
            W_I_IND = onlyInd(m2W) # Indels for immortal (all genes)

            for gene in genes:
                W_M_SNP_G = onlyGene(W_M_SNP, gene)
                W_M_IND_G = onlyGene(W_M_IND, gene)
                W_I_SNP_G = onlyGene(W_I_SNP, gene)
                W_I_IND_G = onlyGene(W_I_IND, gene)
                
                P_ = only(only(only(P, "gn", gene), "m1", m1), "m2", m2) # Protein filtered to specified gene
                assert(len(P_) == 0 or len(P_) == 1)

                P_M_A = "-" if len(P_) == 0 else P_[0]["a1"] # Average mean for mortal samples
                P_I_A = "-" if len(P_) == 0 else P_[0]["a2"] # Average mean for immortal samples

                #
                # Checks mutations and if no mutations add a "-" sign (pathway not activated due to mutations). Otherwise it's "+".
                #
                
                WGS_M_S = gene + "-" if len(W_M_SNP_G) == 0 and len(W_M_IND_G) == 0 else gene + "+" # Gene pathway for mortal
                WGS_I_S = gene + "-" if len(W_I_SNP_G) == 0 and len(W_I_IND_G) == 0 else gene + "+" # Gene pathway for immortal
            
                w.write(f.format(name, m1, m2, gene, len(W_M_SNP_G), len(W_M_IND_G), len(W_I_SNP_G), len(W_I_IND_G), WGS_M_S, WGS_I_S, \
                                 P_M_A, P_I_A))
    w.close()

def parseP(file):
    with open(file, "r") as r:
        for l in r:
            toks = l.strip().split(';')        
            if toks[0] == "Mortal":
                continue
        
            m1 = toks[0] # Mortal
            m2 = toks[1] # Immortal
            gn = toks[3]
            lf = toks[6] # LogFC
            a1 = toks[4] # Mortal mean
            a2 = toks[5] # Immortal mean
            
            if m1 == "IIICF_P7":
                m1 = "IIICF"
            if m2 == "IIICF_P7":
                m2 = "IIICF"
            
            yield { "m1":m1, "m2":m2, "gn":gn, "lf":lf, "a1":a1, "a2":a2 }
            
def parseW(file):
    with open(file, "r") as r:
        for l in r:
            toks = l.strip().split(';')

            # Eg: IIICF-T_B3
            name = file2Samp(toks[1])

            # Translated gene name (from Ensembl)
            gn = ens2Name(toks[19])

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

            # Mutation type
            ty = "SNP" if len(ref) == 1 and len(alt) == 1 else "Ind"

            yield { "name":name, "chr":chr, "pos":pos, "ref":ref, "alt":alt, "imp":imp, "gn":gn, "type":ty }

if sys.argv[1] == "W":
    save("AN1/AN1_W.pickle", list(parseW("AN1/FILTERED_WGS.csv")))
elif sys.argv[1] == "P":
    save("AN1/AN1_P.pickle", list(parseP("SWATH2/SWATH2_results.tsv")))
elif sys.argv[1] == "A":
    analyze(load("AN1/AN1_W.pickle"), load("AN1/AN1_P.pickle"))
elif sys.argv[1] == "F":
    pass
