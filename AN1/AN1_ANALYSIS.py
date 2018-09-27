#
# Python script for analyzing resulting after filtering
#
#     python3 AN1/AN1_ANALYSIS.py F "/home/twong/das/Outbox/Cancer Research Unit/RReddell/WGS/germline.csv" > AN1/AN1_W_FILTERED.csv
#     python3 AN1/AN1_ANALYSIS.py W
#     python3 AN1/AN1_ANALYSIS.py P
#     python3 AN1/AN1_ANALYSIS.py A
#
# Generate a user-friendly TSV output file for the results. Further statisticsl analysis can be done on the output file.
#

import sys
import pickle

# Genes interested
genes = [ "ENSG00000204209", "ENSG00000085224", "ENSG00000164362", "ENSG00000141510", "DAXX", "ATRX", "TERT", "TP53" ]

# DAXX, ATRX, TERT and TP53
rs = [ "chr6:33318558", "chrX:77504878", "chr5:1253147", "chr17:7661779" ]

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
    f = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"

    w = open("AN1/AN1_RESULTS.tsv", "w")
    w.write(f.format("Name", "Mortal", "Immortal", "Gene/Protein", \
                     "WGS_Mortal_SNP", "WGS_Mortal_Ind", "WGS_Immortal_SNP", "WGS_Immortal_Ind", "WGS_Mortal_S", "WGS_Immortal_S"))

    with open("AN1/AN1_CONTRASTS.csv") as r:
        for l in r:
            if "Mortal" in l and "Immortal" in l:
                continue
            toks = l.strip().split(',')

            m1  = toks[0] # Mortal
            m2  = toks[1] # Immortal
            m1W = WGSByName(W, m1)
            m2W = WGSByName(W, m2)
            
            #assert(len(m1W) > 0)
            #assert(len(m2W) > 0)
            
            # Name of the contrast
            name = m1 + "_" + m2

            W_M_SNP = onlySNP(m1W) if len(m1W) > 0 else "-" # SNPs for mortal (all genes)
            W_M_IND = onlyInd(m1W) if len(m1W) > 0 else "-" # Indels for mortal (all genes)
            W_I_SNP = onlySNP(m2W) if len(m2W) > 0 else "-" # SNPs for immortal (all genes)
            W_I_IND = onlyInd(m2W) if len(m2W) > 0 else "-" # Indels for immortal (all genes)

            for gene in genes:
                W_M_SNP_G = onlyGene(W_M_SNP, gene) if W_M_SNP != "-" else "-" 
                W_M_IND_G = onlyGene(W_M_IND, gene) if W_M_IND != "-" else "-" 
                W_I_SNP_G = onlyGene(W_I_SNP, gene) if W_I_SNP != "-" else "-" 
                W_I_IND_G = onlyGene(W_I_IND, gene) if W_I_IND != "-" else "-" 
                
                #
                # Checks mutations and if no mutations add a "-" sign (pathway not activated due to mutations). Otherwise it's "+".
                #
                
                if W_M_SNP_G == "-" or W_M_IND_G == "-":
                    WGS_M_S = "-"
                else:
                    WGS_M_S = gene + "-" if len(W_M_SNP_G) == 0 and len(W_M_IND_G) == 0 else gene + "+" # Gene pathway for mortal
                    
                if W_I_SNP_G == "-" or W_I_IND_G == "-":
                    WGS_I_S = "-"
                else:                
                    WGS_I_S = gene + "-" if len(W_I_SNP_G) == 0 and len(W_I_IND_G) == 0 else gene + "+" # Gene pathway for immortal
            
                w.write(f.format(name, m1, m2, gene, len(W_M_SNP_G), len(W_M_IND_G), len(W_I_SNP_G), len(W_I_IND_G), WGS_M_S, WGS_I_S))
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
        
            # Predict the effect of coding variants on protein function (e.g. deleterious and deleterious_low_confidence)
            sift = toks[48]
            
            # Possible impact of amino acid substitutions on the stability and function of human proteins
            phen = toks[49]
            
            # Impact (LOW, HIGH, MODERATE, MODIFIER)
            imp = toks[17]
            
            # Consequence (e.g. TF_binding_site_variant, 3_prime_UTR_variant, 5_prime_UTR_variant and regulatory_region_variant)
            cons = toks[16]
            
            # Shortest ditance from variant to transcript
            dist = toks[33]
            
            assert(isMod(imp) or isHigh(imp) or isLow(imp))

            chr = toks[11] # Chromosome
            pos = toks[12] # Position
            ref = toks[13] # Reference
            alt = toks[14] # Allele

            # Mutation type
            ty = "SNP" if len(ref) == 1 and len(alt) == 1 else "Ind"

            # How to label this variant? (either inside a gene or a promoter)
            lab = "Normal" if gn in genes else "Promoter"

            yield { "name":name, "lab":lab, "chr":chr, "pos":pos, "ref":ref, "alt":alt, "imp":imp, "gn":gn, "type":ty, "sift":sift, "phen":phen }

if sys.argv[1] == "W":
    save("AN1/AN1_W.pickle", list(parseW("AN1/AN1_W_FILTERED.csv")))
elif sys.argv[1] == "P":
    save("AN1/AN1_P.pickle", list(parseP("SWATH2/SWATH2_results.tsv")))
elif sys.argv[1] == "A":
    analyze(load("AN1/AN1_W.pickle"), load("AN1/AN1_P.pickle"))
elif sys.argv[1] == "F":
    with open(sys.argv[2], "r") as r:
        chrs = [ { "c":i.split(":")[0], "p1":int(i.split(":")[1])-5000, "p2":int(i.split(":")[1]) } for i in rs]
        for line in r:
            toks = line.split(";")
            
            # Gene name
            gn = toks[19]
            
            # Chromosome
            c = toks[11]
            
            # Position
            p = int(toks[12])

            if any(x == gn for x in genes):
                print(line, end='')
            else:
                if any((c == x["c"] and p >= x["p1"] and p <= x["p2"]) for x in chrs):                
                    print(line, end='')                    
