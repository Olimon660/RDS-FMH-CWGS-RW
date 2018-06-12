import sys

with open(sys.argv[1]) as r:
    for line in r:
        toks = line.strip().split('\t') 
        assert(len(toks) == 2)
        normal = toks[0]
        tumor  = toks[1]               
        cmd = 'qsub -v normalID=' + normal + ',tumorID=' + tumor + ' MUTECT_with_control.pbs'
        print(cmd)