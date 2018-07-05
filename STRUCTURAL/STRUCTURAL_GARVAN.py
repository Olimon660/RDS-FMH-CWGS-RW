#
# python STRUCTURAL/STRUCTURAL_GARVAN.py STRUCTURAL/GARVAN.txt
#

import os
import sys

with open(sys.argv[1]) as r:
    for line in r:
        toks = line.strip().split('\t') 
        assert(len(toks) == 2)
        
        tumor  = toks[0]
        normal = toks[1]
        
        cmd = 'qsub -m beas -pe smp 2 -l mem_requested=16G,tmp_requested=16G -q short.q -V -cwd  -j y -S /bin/sh STRUCTURAL/STRUCTURAL_GARVAN.pbs ' + normal + ' ' + tumor
        print(cmd)
        os.system(cmd)
