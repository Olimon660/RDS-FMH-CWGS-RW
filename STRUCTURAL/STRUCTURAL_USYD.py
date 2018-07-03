#
# python STRUCTURAL/STRUCTURAL_USYD.py samples.txt
#

import os
import sys

with open(sys.argv[1]) as r:
    for line in r:
        sample = line.strip()
        
        cmd = 'qsub -v sample=' + sample + ' STRUCTURAL/STRUCTURAL_SOMATIC.pbs'
        print(cmd)
        os.system(cmd)
