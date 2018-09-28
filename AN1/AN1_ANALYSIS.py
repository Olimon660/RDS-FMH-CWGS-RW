#
# Python script for analyzing resulting after filtering
#
#     python3 AN1/AN1_ANALYSIS.py F "/home/twong/das/Outbox/Cancer Research Unit/RReddell/WGS/germline.csv" > AN1/AN1_W_FILTERED.csv
#     python3 AN1/AN1_ANALYSIS.py W
#     python3 AN1/AN1_ANALYSIS.py P
#     python3 AN1/AN1_ANALYSIS.py T
#     python3 AN1/AN1_ANALYSIS.py A
#     python3 AN1/AN1_ANALYSIS.py R
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

def grep(x, key, v):
    return [i for i in x if v in i[key]]

def only(x, key, v):
    return [i for i in x if i[key] == v]

def analyze(W, P):
    w1 = open("AN1/AN1_SUMMARY.tsv", "w")
    
    #
    # - Name
    # - Mortal namae
    # - Immortal name
    # - Feature (gene or protein)
    #
    # - Number of SNPS in WGS mortal
    # - Number of indels in WGS mortal
    # - Number of consequence in WGS mortal (upstream_gene_variant)
    # - Number of consequence in WGS mortal (regulatory_region_variant)
    # - Number of consequence in WGS mortal (TF_binding_site_variant)
    # - Number of consequence in WGS mortal (5_prime_UTR_variant&NMD_transcript_variant or 5_prime_UTR_variant)
    # - Number of impacts in WGS mortal (HIGH)
    # - Number of impacts in WGS mortal (MODERATE)
    # - Number of impacts in WGS mortal (LOW)
    #
    # - Number of SNPS in WGS immortal
    # - Number of indels in WGS immortal
    # - Number of consequence in WGS immortal (upstream_gene_variant)
    # - Number of consequence in WGS immortal (regulatory_region_variant)
    # - Number of consequence in WGS immortal (TF_binding_site_variant)
    # - Number of consequence in WGS immortal (5_prime_UTR_variant&NMD_transcript_variant or 5_prime_UTR_variant)
    # - Number of impacts in WGS immortal (HIGH)
    # - Number of impacts in WGS immortal (MODERATE)
    # - Number of impacts in WGS immortal (LOW)
    #
    
    # Format string
    f = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
    
    w1.write(f.format("Name", "Mortal", "Immortal", "Feature", \
                      "WM_SNP", "WM_Ind", "WM_ConUp", "WM_ConReg", "WM_ConTF", "WM_Con5Prime", "WM_ConHigh", "WM_ConMod", "WM_ConLow",
                      "WI_SNP", "WI_Ind", "WI_ConUp", "WI_ConReg", "WI_ConTF", "WI_Con5Prime", "WI_ConHigh", "WI_ConMod", "WI_ConMod" ))

    with open("AN1/AN1_CONTRASTS.csv") as r:
        for l in r:
            if "Mortal" in l and "Immortal" in l:
                continue
            toks = l.strip().split(',')

            m1 = toks[0] # Mortal
            m2 = toks[1] # Immortal

            # Construct a block of text for mortal/immortal for a gene
            def block(m, gn):
                x = only(only(W, "name", m), "gn", gn)
                #assert(len(x) > 0) Wait until new annotation is done

                snp = only(x, "type", "snp") # SNPs for mortal
                ind = only(x, "type", "ind") # Indels for mortal
                ups = only(x, "con",  "upstream_gene_variant")
                reg = only(x, "con",  "regulatory_region_variant")
                tfs = only(x, "con",  "TF_binding_site_variant")
                fiv = grep(x, "con",  "5_prime_UTR_variant")
                hig = only(x, "imp",  "HIGH")
                med = only(x, "imp",  "MODERATE")
                low = only(x, "imp",  "LOW") + only(x, "imp", "MODIFIER")
                
                return (str(len(snp)) + "\t" + str(len(ind)) + "\t" + str(len(ups)) + "\t" + str(len(reg)) + "\t" + str(len(tfs)) + "\t" + str(len(fiv)) + \
                        str(len(hig)) + "\t" + str(len(med)) + "\t" + str(len(low)))
            
            # Write a gene for each contrast 
            for gn in genes:
                w1.write((m1 + "_" + m2) + "\t" + m1 + "\t" + m2 + "\t" + gn + "\t" + block(m1, gn) + "\t" + block(m2, gn))
    w1.close()

def parseP(file):
    with open(file, "r") as r:
        for l in r:
            toks = l.strip().split(';')        
            if toks[0] == "Mortal":
                continue
        
            m1 = toks[0] # Mortal
            m2 = toks[1] # Immortal
            gn = toks[3] # Ensembl gene
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
            
            ft = toks[20] # Feature type
            ff = toks[21] # Feature
            
            sift = toks[48] # Predict the effect of coding variants on protein function (e.g. deleterious and deleterious_low_confidence)
            phen = toks[49] # Possible impact of amino acid substitutions on the stability and function of human proteins            

            imp = toks[17] # Impacts
            con = toks[16] # Consequence
            
            # Shortest ditance from variant to transcript
            dist = toks[33]

            chr = toks[11] # Chromosome
            pos = toks[12] # Position
            ref = toks[13] # Reference
            alt = toks[14] # Allele

            # Mutation type
            ty = "snp" if len(ref) == 1 and len(alt) == 1 else "ind"

            # How to label this variant? (either inside a gene or a promoter)
            lab = "Normal" if gn in genes else "Promoter"

            yield { "name":name, "lab":lab, "chr":chr, "pos":pos, "ref":ref, "alt":alt, "con":con, "imp":imp, "gn":gn, "type":ty, "sift":sift, "phen":phen }

if sys.argv[1] == "W":
    save("AN1/AN1_W.pkl", list(parseW("AN1/AN1_W_FILTERED.csv")))
elif sys.argv[1] == "S":
    pass
elif sys.argv[1] == "A":
    analyze(load("AN1/AN1_W.pkl"), None)
elif sys.argv[1] == "F":
    with open(sys.argv[2], "r") as r:
        # Make sure we capture all upstream variants
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
            elif c == "chrE" or "chrC" or "chrV": # Any of the special genome
                print(line, end='')
            else:
                if any((c == x["c"] and p >= x["p1"] and p <= x["p2"]) for x in chrs):                
                    print(line, end='')                    
