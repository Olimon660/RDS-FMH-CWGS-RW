#
# python3 Somatic/SOMATIC_SQL.py 7 somaticAnnotated.bin
#

import os
import sys
import sqlite3

def parseTxt(file, name, c):
    with open(file, 'r') as f:
        for line in f:
            if len(line) > 0 and line[0] == '#':
                continue
            toks = line.split('\t')
            chrom = toks[1].split(':')[0]
            start = toks[1].split(':')[1]
            ref   = toks[0].split('_')[2].split('/')[0]
            alt   = toks[0].split('_')[2].split('/')[1]
            
            c.execute('INSERT INTO Variants(Name, Chrom, Start, Ref, Alt, Gene, Feature, Type, Info)' \
                      'VALUES("{name}", "{chrom}", "{start}", "{ref}", "{alt}", "{gene}", "{feature}", "{type}", "{info}")'.format(\
                       name=name, chrom=chrom, start=start, ref=ref, alt=alt, gene=toks[3],\
                       feature=toks[4], type=toks[5], info=toks[6]))

SQL = """CREATE TABLE Variants (
            Name  varchar(255),
            Chrom varchar(255),
            Start varchar(255),
            Ref   varchar(255),
            Alt   varchar(255),
            Gene  varchar(255),
            Feature varchar(255),
            Type    varchar(255),
            Info    varchar(255)
       );"""

os.system('rm -f ' + sys.argv[2])
conn = sqlite3.connect(sys.argv[2])
c = conn.cursor()
c.execute(SQL)

for (dirpath, dirs, files) in os.walk(sys.argv[1]):
    for file in files:
        if 'ANNOTATED_REMOVED' in file and file.endswith('.vcf'):
            txt  = sys.argv[1] + os.sep + file
            html = txt + '_summary.html'
            warn = txt + '_warnings.txt'

            assert(os.path.isfile(txt))
            assert(os.path.isfile(html))
            assert(os.path.isfile(warn))                        
            
            txt = parseTxt(txt, file, c)            
            conn.commit()