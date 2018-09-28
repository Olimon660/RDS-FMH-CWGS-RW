#
# python3 GERMLINE/GERMLINE_CSV.py OTHERS/ANNOTATION.sql <Input File> <Output CSV>
#

import os
import sys
import vcf
import csv
import json
import sqlite3

def insert(c, conn, vs, cols):
    assert(len(vs) > 0)
    VAL = ''
    for i in range(len(cols)):
        if i != len(cols) - 1:
            VAL = VAL + '\'%s' + '\','
        else:
            VAL = VAL + '\'%s' + '\''
            
    for i in vs:
        SQL = 'INSERT INTO Annotation VALUES('
        for j in range(len(cols)):
            x = i[cols[j]]
            if x is None or str(x) == '':
                SQL = SQL + ' NULL' + ','
            else:            
                SQL = SQL + ' "' + str(x) + '",'
        SQL += (')')
        SQL = SQL.replace(",)", ')').replace(' "', '"')
        c.execute(SQL)

    conn.commit()
    
def parseSQL(file):
    x = []
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if 'VARCHAR' in line:
                toks = [i for i in line.replace('\t', ' ').split(' ') if i != '']
                assert('VARCHAR' in toks[1])
                x.append(toks[0])
    return x

def parseHead(file):
    with open(file) as r:
        for line in r:
            line = line.strip()
            if line[0] != '#':
                break
            elif 'CSQ' in line:
                toks = line.split('. Format: ')
                return toks[1].replace('">', '').split('|')

def parseVCF(file, csv, cols):
    h = parseHead(file)
    assert(len(h) > 0)

    r = vcf.Reader(open(file, 'r'))
    n = 0

    for i in r:
        key = str(i.CHROM) + '_' + str(i.POS) + '_' + str(i.REF) + '/' + str(i.ALT[0])
        n = n + 1
        if (n % 10000) == 0:
            print(n)

        DP = i.INFO['DP'] # Total read depth

        def AD1():
            if 'AD' in i.FORMAT:
                return i.samples[0]['AD']
            else:
                return None
                
        def AD2():
            if 'AD' in i.FORMAT and len(i.samples) > 1:
                return i.samples[1]['AD']
            else:
                return None
            
        AD_1 = AD1()
        AD_2 = AD2()
        
        def AF1():
            if 'AF' in i.INFO:
                return i.INFO['AF']
            else:
                return i.samples[0]['AF'] # Allele frequency for sample 1

        def AF2():
            return None if len(i.samples) == 1 else i.samples[1]['AF'] # Allele frequency for sample 2

        def GT1():
            if 'GT' in i.INFO:
                return i.INFO['GT']
            else:
                return i.samples[0]['GT']

        def GT2():
            return None if len(i.samples) == 1 else i.samples[1]['GT']
        
        def QUAL1():
            if 'NLOD' in i.INFO:
                return str(i.INFO['NLOD'][0]) + ':' + str(i.INFO['N_ART_LOD'][0]) # Somatic analysis
            else:
                return i.QUAL                

        def QUAL2():
            return None if len(i.samples) == 1 else str(i.INFO['TLOD'][0])            
            
        AF_1 = AF1()
        AF_2 = AF2()
        GT_1 = GT1()
        GT_2 = GT2()        
        QL_1 = QUAL1() # Quality for sample 1 (could be multi-string)
        QL_2 = QUAL2() # Quality for sample 2 (could be multi-string)

        for j in range(len(i.INFO['CSQ'])):
            CSQ = i.INFO['CSQ'][j]
            toks = CSQ.split('|')
            assert(len(h) == len(toks))
            x = [ key, file, DP, AD_1, AD_2, AF_1, AF_2, GT_1, GT_2, QL_1, QL_2, str(i.CHROM), str(i.POS), i.REF, i.ALT[0]]
            for k in range(len(toks)):
                x.append(toks[k])
            assert(len(x) == len(cols))
            csv.writerow(x)

cols = parseSQL(sys.argv[1])
csv  = csv.writer(open(sys.argv[3], 'w'), delimiter=';', quoting=csv.QUOTE_MINIMAL)

# Directory of annotated and filtered VCF files?
for (dirpath, dirs, files) in os.walk(sys.argv[2]):
    for file in files:
        if 'ANNOTATED_REMOVED' in file and file.endswith('.vcf'):
            print(file)
            parseVCF(sys.argv[2] + os.sep + file, csv, cols)
