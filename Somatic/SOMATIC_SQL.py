#
# python3 Somatic/SOMATIC_SQL.py 7 database.bin
#

import os
import sys
import sqlite3

def insert(c, x):
    SQL = 'INSERT INTO Annotation VALUES("%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s")'
    SQL = SQL % (x['Key'], x['File'], x['Chrom'], x['Start'], x['REF'], x['ALT'], x['Type'], x['Gene'], \
                 x['Feature'], x['FeatureType'], x['Consequence'], x['cDNAPosition'], x['CDSPosition'], \
                 x['ProteinPosition'], x['AminoAcids'], x['Codons'], x['ExistingVariation'], x['Extra'])
    c.execute(SQL)

def parse(file, name):
    with open(file, 'r') as f:
        for line in f:
            if len(line) > 0 and line[0] == '#':
                continue
            toks = line.split('\t')
            
            def NULL(x):
                return None if x == '-' else x
            
            REF = toks[0].split('_')[2].split('/')[0]
            ALT = toks[0].split('_')[2].split('/')[1]
            
            for i in toks[13].split(';'):
                key = i.split('=')[0]
                val = i.split('=')[1]

            info  = { 'Key'  : toks[0], \
                      'File' : file,    \
                      'Chrom': toks[1].split(':')[0], \
                      'Start': toks[1].split(':')[1], \
                      'REF'  : REF, \
                      'ALT'  : ALT, \
                      'Type' : 'SNP' if len(REF) == len(ALT) else 'Indel', \
                      'Gene' :             toks[3],  \
                      'Feature':           toks[4],  \
                      'FeatureType':       toks[5],  \
                      'Consequence':       toks[6],  \
                      'cDNAPosition':      NULL(toks[7]),  \
                      'CDSPosition':       NULL(toks[8]),  \
                      'ProteinPosition':   NULL(toks[9]),  \
                      'AminoAcids':        NULL(toks[10]), \
                      'Codons':            NULL(toks[11]), \
                      'ExistingVariation': NULL(toks[12]), \
                      'Extra':             NULL(toks[13])  \
            }
            
            insert(c, info)

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

            txt = parse(txt, file)            
            conn.commit()