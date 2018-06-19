import os
import sys

with open(sys.argv[1]) as r:
    for line in r:
        toks = line.strip().split('\t') 
        assert(len(toks) == 2)
        
        #
        # The file follows Erdahl's Excel. IMMORTAL is on first column and is tumor, while MORTAL is normal on the
        # second column.
        #
        
        tumor  = toks[0]
        normal = toks[1]
        
        cmd = 'qsub -v normalID=' + normal + ',tumorID=' + tumor + ' STAGE2/gridss.pbs'
        print(cmd)
        #cmd = 'qsub -v normalID=' + normal + ',tumorID=' + tumor + ' STAGE2/TEST.pbs'        

        #os.system(cmd)
