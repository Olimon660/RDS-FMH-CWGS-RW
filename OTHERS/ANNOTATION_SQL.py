#
# python3 OTHERS/ANNOTATION_SQL.py 7 database.bin OTHERS/ANNOTATION.sql 
#

import os
import sys
import vcf
import json
import sqlite3

def insert(c, vs, cols):
    VAL = ''
    for i in range(len(cols)):
        if i != len(cols) - 1:
            VAL = VAL + '\'%s' + '\','
        else:
            VAL = VAL + '\'%s' + '\''
    for i in vs:
        SQL = 'INSERT INTO Annotation VALUES('
        for j in range(len(cols)):
            SQL = SQL + ' "' + str(i[cols[j]]) + '",'
        SQL += (')')
        SQL = SQL.replace(",)", ')').replace(' "', '"')
        c.execute(SQL)

def parseHead(file):
    with open(file) as r:
        for line in r:
            line = line.strip()
            if line[0] != '#':
                break
            elif 'CSQ' in line:
                toks = line.split('. Format: ')
                return toks[1].replace('">', '').split('|')

def parseVCF(file, cols):
    h = parseHead(file)
    assert(len(h) > 0)
    x = []
    r = vcf.Reader(open(file, 'r'))

    for i in r:
        key = str(i.CHROM) + '_' + str(i.POS) + '_' + str(i.REF) + '/' + str(i.ALT[0])

        DP = i.INFO['DP'] # Total read depth
        
        AD_1 = i.samples[0]['AD'] # Allele read depth for sample 1
        AD_2 = i.samples[1]['AD'] if len(i.samples) > 1 else None        
        AF_1 = i.samples[0]['AF'] # Allele frequency for sample 1
        AF_2 = i.samples[1]['AF'] if len(i.samples) > 1 else None
        GT_1 = i.samples[0]['GT'] # Genotype for sample 1
        GT_2 = i.samples[1]['GT'] if len(i.samples) > 1 else None        

        def QUAL1(x):
            return str(x['NLOD'][0]) + ':' + str(x['N_ART_LOD'][0])

        def QUAL2(x):
            return str(x['TLOD'][0])
        
        QL_1 = QUAL1(i.INFO) # Quality for sample 1 (could be multi-string)
        QL_2 = QUAL2(i.INFO) # Quality for sample 2 (could be multi-string)
        
        CSQ = i.INFO['CSQ'][0]
        toks = CSQ.split('|')
        assert(len(h) == len(toks))
        
        v = { 'File':file,          \
              'Chrom':str(i.CHROM), \
              'POS':str(i.POS),     \
              'REF':i.REF,          \
              'ALT':i.ALT[0],       \
              'Key':key, 'DP':DP, 'AD_1':AD_1, 'AD_2':AD_2, 'AF_1':AF_1, 'AF_2':AF_2, 'GT_1':GT_1, 'GT_2':GT_2, 'QL_1':QL_1, 'QL_2':QL_2 }

        for j in range(len(toks)):
            v[h[j]] = toks[j]
        x.append(v)

    return x
    
def parseSQL(file):
    x = []
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if 'VARCHAR' in line:
                toks = [i for i in line.replace('\t', ' ').split(' ') if i != '']
                assert(len(toks) == 2 or len(toks) == 4)
                x.append(toks[0])
    return x

os.system('rm -f ' + sys.argv[2])
conn = sqlite3.connect(sys.argv[2])
c = conn.cursor()
with open(sys.argv[3], 'r') as f:
    c.executescript(f.read())
cols = parseSQL(sys.argv[3])

for (dirpath, dirs, files) in os.walk(sys.argv[1]):    
    for file in files:
        if 'ANNOTATED_REMOVED' in file and file.endswith('.vcf'):
            print(file)
            insert(c, parseVCF(sys.argv[1] + os.sep + file, cols), cols)
            conn.commit()