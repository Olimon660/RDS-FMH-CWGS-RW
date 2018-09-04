#
# python3 SWATH2/SWATH2_SAMPLE.py /tmp/A.txt /tmp/B.txt
#

import sys

x = []
with open(sys.argv[1], "r") as r:
    for line in r:
        line = line.strip().replace(" ", "_")
        assert(len(line.split(" ")) == 1)
        assert(len(line.split("\t")) == 1)
        if "Cell_Study_Sample" in line and not "Cell_Study_Sample_" in line:
            line = line.replace("Cell_Study_Sample", "Cell_Study_Sample_")
        line = line.replace("Sample0", "Sample_0")
        line = line.replace("Sample1", "Sample_1")
        line = line.replace("Sample2", "Sample_2")
        line = line.replace("Sample3", "Sample_3")
        line = line.replace("Sample4", "Sample_4")
        line = line.replace("Sample5", "Sample_5")
        line = line.replace("Sample6", "Sample_6")
        line = line.replace("Sample7", "Sample_7")
        line = line.replace("Sample8", "Sample_8")
        line = line.replace("Sample9", "Sample_9")
        line = line.replace("_Sample.", "_Sample_")
        line = line.replace("..Cell.Study_Sample.wiff..sample.1..", "")        
        for i in range(0,10):
            line = line.replace("..Cell.Study_Sample" + str(i) + ".wiff..sample.1..", "")

        if "_" in line:
            if ".wiff" in line:
                line = line.split(".wiff")[0]            
            line = (line.split("_"))[-1]

        assert(line.isdigit())
        x.append(int(line))

assert(len(x) == len(set(x)))

with open(sys.argv[2], "w") as w:
    for i in x:
        w.write(str(i) + "\n")
