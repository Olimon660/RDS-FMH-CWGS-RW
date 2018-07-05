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
        
        cmd = 'STRUCTURAL/STRUCTURAL_GARVAN.pbs ' + normal + ' ' + tumor
        cmd = 'qsub -q short.q -m beas -pe smp 2 -l mem_requested=16G,tmp_requested=16G -V -cwd -j y -S ' + cmd
        
        qsub -q short.q -pe smp 8 -l mem_requested=16G,tmp_requested=60G COMMAND
        
        print(cmd)
        os.system(cmd)
