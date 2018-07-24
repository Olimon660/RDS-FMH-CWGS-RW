#
# python3 STRUCTURAL/STRUCTURAL_MERGE.py 8/STRUCTURAL_SQL.csv
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

            t = {}
            for i in c:
                t[i] = toks[c[i]]
            x1.append(t)
            
            key1 = toks[c['Sample']] + '_' + toks[c['Name']]
            key2 = toks[c['Sample']] + '_' + toks[c['Partner']]

            x2[key1] = t
            x2[key2] = t 

    return (x1, x2)

(b1, b2) = readTSV('8/STRUCTURAL_SUBTRACT.tsv')
csv = open(sys.argv[1], 'w')

o1 = [ 'Type', 'File', 'Sample', 'Resolu', 'SVType', 'Chr1', 'Start1', 'End1', 'Chr2', 'Start2', 'End2', 'Ref1', 'Alt1',   \
       'Ref2', 'Alt2', 'SVLen', 'RP1', 'RP2', 'RPQ1', 'RPQ2', 'SR1', 'SR2', 'SRQ1', 'SRQ2', 'REFC1', \
       'REFC2', 'REFPAIR1', 'REFPAIR2', 'HOMLEN1', 'HOMSEQ1', 'HOMLEN2', 'HOMSEQ2' ]
o2 = [ 'Allele', 'Annotation', 'AnnotationImpact', 'GeneName', 'GeneID', 'FeatureType', 'FeatureID', \
       'TranscriptBioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA', 'CDS', 'AA', 'Distance', 'Status' ]

for file in os.listdir("8"):
    if file.startswith('SUBTRACTED') and file.endswith(".vcf"):
        print(file)
        r = vcf.Reader(open(os.path.join("8", file), 'r'))        
        for i in r:
            samp = file.replace('SUBTRACTED_ANNOTATED_SOMATIC_GRIDSS_', '').replace('.vcf', '')
            name = i.ID
            key1 = samp + '_' + name
            
            # Has the key been filtered out by normal?
            if key1 in b2:
                b = b2[key1]
                s1 = ''
                for o in o1:
                    if s1 == '':
                        s1 = s1 + b[o]
                    else:
                        s1 = (s1 + ';' + b[o])                

                if not 'ANN' in i.INFO:
                    s2 = ''
                    for o in o2:
                        s2 = (s2 + ';' + '-')
                    csv.write(s1 + s2 + '\n')                        
                else:
                    for a in i.INFO['ANN']:
                        ann = a.split('|')
                        assert(len(ann) == 16)
                        s2 = ''
                        for k in range(0, len(o2)):
                            s2 = (s2 + ';' + ann[k])
                        assert(len(s1.split(';')) == len(o1))
                        assert(len(s2.split(';')) == len(o2) + 1)
                        csv.write(s1 + s2 + '\n')
            else:
                raise Exception('Unknown key: ' + key1)
