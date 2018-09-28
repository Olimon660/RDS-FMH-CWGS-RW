#
# Python script for analyzing resulting after filtering
#
#     python3 AN1/AN1_ANALYSIS.py F "/home/twong/das/Outbox/Cancer Research Unit/RReddell/WGS/germline.csv" > AN1/AN1_W_FILTERED.csv
#     python3 AN1/AN1_ANALYSIS.py G
#     python3 AN1/AN1_ANALYSIS.py V
#     python3 AN1/AN1_ANALYSIS.py A
#     python3 AN1/AN1_ANALYSIS.py R
#
# Generate a user-friendly TSV output file for the results. Further statisticsl analysis can be done on the output file.
#

import sys
import pickle

# Genes interested
gns = [ "ENSG00000204209", "ENSG00000085224", "ENSG00000164362", "ENSG00000141510", "DAXX", "ATRX", "TERT", "TP53" ]

# DAXX, ATRX, TERT and TP53
rs = [ "chr6:33318558", "chrX:77504878", "chr5:1253147", "chr17:7661779" ]

# Eg: 6/ANNOTATED_REMOVED_FILTERED_INDEL_NORM_DECOM_GATK_IIICF-T_B3.vcf
def file2Samp(x):
    assert("GATK" in x)
    return x.split("GATK_")[1].replace(".vcf", "").replace("-", "_")

def ens2Name(x):
    if x == "ENSG00000204209":
        return "DAXX"
    elif x == "ENSG00000085224":
        return "ATRX"
    elif x == "ENSG00000164362":
        return "TERT"
    elif x == "ENSG00000141510":
        return "TP53"
    else:
        return "-"        
        #raise Exception("Unknown: " + str(x))

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

def rate(file):
    pass

def analyze(W, V):
    print("Writing AN1/AN1_SUMMARY.tsv")
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
    # - Number of mutations in SV mortal
    # - Number of mutations in SV immportal
    # - Number of impacts in WGS mortal (HIGH)
    # - Number of impacts in WGS mortal (MODIFIER)
    # - Number of impacts in WGS mortal (LOW)
    # - Number of impacts in WGS immortal (HIGH)
    # - Number of impacts in WGS immortal (MODIFIER)
    # - Number of impacts in WGS immortal (LOW)
    #
    
    w1.write("Name\tMortal\tImmortal\tFeature\t"
             "GM_SNP\tGM_Ind\tGM_ConUp\tGM_ConReg\tGM_ConTF\tGM_Con5Prime\tGM_ImpHigh\tGM_ImpMod\tGM_ImpLow\t"
             "GI_SNP\tGI_Ind\tGI_ConUp\tGI_ConReg\tGI_ConTF\tGI_Con5Prime\tGI_ImpHigh\tGI_ImpMod\tGI_ImpLow\t"
             "VM_MUT\tVI_MUT\tVM_ImpHigh\tVM_ImpMod\tVM_ImpLow\tVI_ImpHigh\tVI_ImpMod\tVI_ImpLow\n")

    with open("AN1/AN1_CONTRASTS.csv") as r:
        for l in r:
            if "Mortal" in l and "Immortal" in l:
                continue
            toks = l.strip().split(',')

            m1 = toks[0] # Mortal
            m2 = toks[1] # Immortal

            # Construct a block of text for mortal/immortal for a gene (germline)
            def blockG(m, gn):
                m = "IIICF_T_C3"
                gn = "TERT"
                
                
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
                
                return (str(len(snp)) + "\t" + str(len(ind)) + "\t" + str(len(ups)) + "\t" + str(len(reg)) + "\t" + str(len(tfs)) + "\t" + str(len(fiv)) + "\t" + \
                        str(len(hig)) + "\t" + str(len(med)) + "\t" + str(len(low)))
            
            # Construct a block of text for mortla/immortal for a gene (SV)
            def blockV(m, gn):                
                x = only(V, "g1", gn)
                
                hig = only(x, "imp", "HIGH")
                med = only(x, "imp", "MODIFIER")
                low = only(x, "imp", "LOW")
                
                return (str(len(x)) + "\t" + str(len(hig)) + "\t" + str(len(med)) + "\t" + str(len(low)))
            
            # Write a gene for each contrast 
            for gn in gns:
                if "ENSG" in gn: # Only the actual gene names
                    continue
                
                w1.write((m1 + "_" + m2) + "\t" + m1 + "\t" + m2 + "\t" + gn + "\t" + blockG(m1, gn) + "\t" + blockG(m2, gn) + "\t" + \
                          blockV(m1, gn) + "\t" + blockV(m2, gn) + "\n")
    w1.close()

#
# Parse input data for proteomics
#

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

#
# Parse input data for structural variant
#

def parseV(file):
    with open(file, "r") as r:
        for l in r:
            toks = l.strip().split(';')
            
            c1 = toks[5]
            s1 = toks[6]
            e1 = toks[7]
            c2 = toks[8]
            s2 = toks[9]
            e2 = toks[10]
            g1 = toks[35] # Gene name
            g2 = toks[36] # Gene ID
            di = toks[46] # Distance            
            ft = toks[37] # Feature type
            fi = toks[38] # Feature ID
    
            # Impacts (-, "LOW", "MODIFIER" and "HIGH")
            imp = toks[34]

            # Only precise SVs are needed
            if toks[3] != "Precise":
                continue
        
            if g1 in gns or g2 in gns:
                yield { "name":toks[1], "type":toks[2], "res":toks[3], "styp":toks[4], \
                        "c1":c1, "s1":s1, "e1":e1, "c2":c2, "s2":s2, "e2":e2, \
                        "g1":g1, "g2":g2, "di":di, "ft":ft, "fi":fi, "imp":imp }
            
#
# Parse input data for germline
#

def parseG(file):
    with open(file, "r") as r:
        for l in r:
            toks = l.strip().split(';')

            # Eg: IIICF-T_B3
            name = file2Samp(toks[1])
            
            print(name)

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
            lab = "Normal" if gn in gns else "Promoter"

            yield { "name":name, "lab":lab, "chr":chr, "pos":pos, "ref":ref, "alt":alt, "con":con, "imp":imp, "gn":gn, "type":ty, "sift":sift, "phen":phen }

if sys.argv[1] == "G":
    save("AN1/AN1_G.pkl", list(parseG("AN1/AN1_W_FILTERED.csv"))) # Germline variants
elif sys.argv[1] == "V":
    save("AN1/AN1_V.pkl", list(parseV("7/7.csv"))) # Structural variants
elif sys.argv[1] == "A":
    analyze(load("AN1/AN1_G.pkl"), load("AN1/AN1_V.pkl"))
elif sys.argv[1] == "F":
    with open(sys.argv[2], "r") as r: # Filtering germline variants
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
            
            if any(x == gn for x in gns):
                print(line, end='')
            elif c == "chrE" or "chrC" or "chrV": # Any of the special genome
                print(line, end='')
            else:
                if any((c == x["c"] and p >= x["p1"] and p <= x["p2"]) for x in chrs):                
                    print(line, end='')                    
