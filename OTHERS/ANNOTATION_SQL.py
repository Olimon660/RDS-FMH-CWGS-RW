#
# python3 OTHERS/ANNOTATION_SQL.py OTHERS/ANNOTATION.sql <Input File> <Output Database>
#

import os
import sys
import vcf
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

def parseVCF(file, c, con, cols):
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

        x = { 'File':file,          \
              'Chrom':str(i.CHROM), \
              'POS':str(i.POS),     \
              'REF':i.REF,          \
              'ALT':i.ALT[0],       \
              'Key':key, 'DP':DP, 'AD_1':AD_1, 'AD_2':AD_2, 'AF_1':AF_1, 'AF_2':AF_2, 'GT_1':GT_1, 'GT_2':GT_2, 'QL_1':QL_1, 'QL_2':QL_2 }

        for j in range(len(i.INFO['CSQ'])):
            CSQ = i.INFO['CSQ'][j]
            toks = CSQ.split('|')
            assert(len(h) == len(toks))        
            for k in range(len(toks)):
                x[h[k]] = toks[k]        
            insert(c, con, [x], cols) # Insert for each CSQ combination
    
os.system('rm -f ' + sys.argv[3])
con = sqlite3.connect(sys.argv[3])
c = con.cursor()
with open(sys.argv[1], 'r') as f:
    c.executescript(f.read())

# Columns from schema
cols = parseSQL(sys.argv[1])

if os.path.isdir(sys.argv[2]):
    for (dirpath, dirs, files) in os.walk(sys.argv[2]):
        for file in files:
            if 'ANNOTATED_REMOVED' in file and file.endswith('.vcf'):
                print(file)
                parseVCF(sys.argv[1] + os.sep + file, c, con, cols)
else:
    parseVCF(sys.argv[2], c, con, cols)

