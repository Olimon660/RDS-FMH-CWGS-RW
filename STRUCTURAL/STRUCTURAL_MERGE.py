#
# python3 STRUCTURAL/STRUCTURAL_SQL.py 8/STRUCTURAL_SQL.csv
#

import os
import csv
import sys
import vcf

def readTSV(file):
    x1 = []
    x2 = {}
    with open(file, "r") as f:
        c = {}
        for line in f:
            toks = line.strip().split('\t')
            if 'Chr1' in line:
                for i in range(len(toks)):
                    c[toks[i]] = i
                continue

            t = { 'SVType':  toks[0],  'Chr1':toks[1],  'Start1':toks[2], 'End1':toks[3],   'Chr2':toks[4],  \
                  'Start2':  toks[5],  'End2':toks[6],  'Ref1':  toks[7], 'Alt1':toks[8],   'Name':toks[9],  \
                  'Partner': toks[10], 'Ref2':toks[11], 'Alt2':toks[12],  'SVLen':toks[13], 'File':toks[14], \
                  'Sample':  toks[15] }

            x1.append(t)
            
            key1 = toks[15] + '_' + toks[9]
            key2 = toks[15] + '_' + toks[10]
            
            if key1 not in x2:
                x2[key1] = []
            if key2 not in x2:
                x2[key2] = []
                
            x2[key1].append(x2)
            x2[key2].append(x2)            
            
    return (x1, x2)

# Nothing annotated here
(b1, b2) = readTSV('8/STRUCTURAL_SUBTRACT.tsv')

csv = csv.writer(open(sys.argv[1], 'w'), delimiter=';', quoting=csv.QUOTE_MINIMAL)
csv.write('Resolution\tSVType\tChr1\tStart1\tEnd1\tChr2\tStart2\tEnd2\tRef1\tAlt1\tRef2\tAlt2\tSVLen\tRP1\tRP2\tRPQ1\tRPQ2\tSR1\tSR2\tSRQ1\tSRQ2\tREFC1\tREFC2\tREFPAIR1\tREPAIR2')

for file in os.listdir("8"):
    if file.startswith('SUBTRACTED') and file.endswith(".vcf"):
        r = vcf.Reader(open(os.path.join("8", file), 'r'))        
        for i in r:
            samp = file.replace('SUBTRACTED_ANNOTATED_SOMATIC_GRIDSS_', '').replace('.vcf', '')
            name = i.ID
            key1 = samp + '_' + name
            
            # Has the key been filtered out by normal?
            if key1 in b2:
                ann = i.INFO['ANN'][1].split('|') # Tumor
                
                # Annoation we'll write out for this somatic variant
                ann = { 'Allele':      ann[0],
                        'Annotation':  ann[1],
                        'AnnotationImpact':ann[2],
                        'GeneName':    ann[3],
                        'GeneID':      ann[4],
                        'FeatureType': ann[5],
                        'FeatureID':   ann[6],
                        'TranscriptBioType': ann[7],
                        'Rank':   ann[8],
                        'HGVS.c': ann[9],
                        'HGVS.p': ann[10],
                        'cDNA':   ann[11],
                        'CDS':    ann[12],
                        'AA':     ann[13],
                        'Distance': ann[14] }
                
                # Breakpoints from the R-package
                b = b2[key1]
                
                
                
                sadaasd
                csv.writerow(x)


            #key = sample + '_' + 
            
            
            






