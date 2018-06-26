#
# python3 Somatic/SOMATIC_SQL.py 7 somaticAnnotated.bin
#

import os
import sys
import sqlite3

def insertVariant(c, x):
    c.execute('INSERT INTO Variant VALUES("{Name}", "{File}", "{Chrom}", "{Start}", "{REF}", "{ALT}", "{isSNP}")' \
              .format(Name=x['Name'], File=x['File'], Chrom=x['Chrom'], Start=x['Start'], REF=x['REF'], ALT=x['ALT'], isSNP=x['isSNP']))

def parseTxt(file, name, c):    
    v = {}
    t = {}
    with open(file, 'r') as f:
        for line in f:
            if len(line) > 0 and line[0] == '#':
                continue
            toks  = line.split('\t')
            name  = toks[0]
            chrom = toks[1].split(':')[0]
            start = toks[1].split(':')[1]
            REF   = toks[0].split('_')[2].split('/')[0]
            ALT   = toks[0].split('_')[2].split('/')[1]
            gene  = toks[3]
            
            if not name in v:
                v[name] = { 'Name': name, \
                            'File': file, \
                            'Chrom':chrom, \
                            'Start':start, \
                            'REF':  REF, \
                            'ALT':  ALT, \
                            'Gene': gene, \
                            'isSNP': 'T' if len(REF) == len(ALT) else 'F' }
                            
            if toks[5] == 'Transcript':
                insert
            else:
                raise Exception('Unknown feature type: ' + toks[5])

    [insertVariant(c, v[i]) for i in v]

os.system('rm -f ' + sys.argv[2])
conn = sqlite3.connect(sys.argv[2])
c = conn.cursor()
with open('schema.sql', 'r') as myfile:
    data=myfile.read()   
    c.executescript(data)

for (dirpath, dirs, files) in os.walk(sys.argv[1]):
    for file in files:
        if 'ANNOTATED_REMOVED' in file and file.endswith('.vcf'):
            print(file)
            txt  = sys.argv[1] + os.sep + file
            html = txt + '_summary.html'
            warn = txt + '_warnings.txt'

            assert(os.path.isfile(txt))
            assert(os.path.isfile(html))
            assert(os.path.isfile(warn))                        
            
            txt = parseTxt(txt, file, c)            
            conn.commit()
            
            break