#
# python3 SWATH_CLEANUP.py
#

import os

cns = {}
new = []

with open("/Users/twong/Sources/RDS-FMH-CWGS-RW/SWATH/sample.txt", "r") as r:
    for line in r:
        samp = line.strip()
        
        if not samp in cns:
            cns[samp] = 1
        else:
            cns[samp] = cns[samp] + 1
        
        if samp == "IIICF/a2_A":
            new.append("IIICF/a2_4")            
        elif samp.endswith("_A"):
            new.append(samp.replace("_A", "_1"))
        elif samp.endswith("_B"):
            new.append(samp.replace("_B", "_2"))
        elif samp.endswith("_C"):
            new.append(samp.replace("_C", "_3"))
        else:
            new.append(samp + "_" + str(cns[samp]))

with open("/Users/twong/Sources/RDS-FMH-CWGS-RW/SWATH/newSample.txt", "w") as w:
    for samp in new:
        w.write(samp + "\n")