#
# python3 SWATH/SWATH_PROCESS.py > SWATH/swath.csv
#

import os

with open("SWATH/diffResults.csv") as r:
    for line in r:
        toks = line.strip().split(';')
        prot = toks[0].split("|")
        gene = prot[2].split("_HUMAN")[0]
        toks[0] = prot[0] + "|" + prot[1] + gene
        line = ";".join(toks).replace("_HUMAN", "")        
        print(line)
