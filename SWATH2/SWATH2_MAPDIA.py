#
# Python script for running MapDIA on SWATH analysis
#
#   python3 SWATH2/SWATH2_MAPDIA.py
#

import sys

trks = {} # Indexed by sample names
with open("SWATH2/SWATH2_track.tsv", "r") as r:
    for line in r:
        toks = line.strip().split('\t')
        if toks[0] != "IntenSample":
            if not toks[2] in trks:
                trks[toks[2]] = []            
            trks[toks[2]].append({ \
                "ID":     int(toks[1]), \
                "Sample": toks[2]
            })
assert(len(trks) > 0)

sams = {}
ints = {}
with open("SWATH2/SWATH2_intensity.tsv", "r") as r:
    for line in r:
        toks = line.strip().split('\t')
        if len(ints) == 0:
            for i in range(0, len(toks)):
                sams[i] = int(toks[i])
                ints[int(toks[i])] = []
        else:
            for i in range(0, len(toks)):
                ints[sams[i]].append(toks[i])

with open("SWATH2/SWATH2_tests.tsv", "r") as r:
    for line in r:
        toks = line.strip().split('\t')

        mo = toks[0] # Mortal (control)
        im = toks[1] # Immortal (treatment)        
        print("Mortal: " + mo + ". Immortal: " + im)

        assert(mo in trks)
        assert(im in trks)
        assert(len(trks[mo]) == 3)
        assert(len(trks[im]) == 3)
        
        assert(trks[mo][0]["ID"] in ints)
        assert(trks[im][0]["ID"] in ints)
        assert(trks[mo][1]["ID"] in ints)
        assert(trks[im][1]["ID"] in ints)
        assert(trks[mo][2]["ID"] in ints)
        assert(trks[im][2]["ID"] in ints)

        
        
        
        break
