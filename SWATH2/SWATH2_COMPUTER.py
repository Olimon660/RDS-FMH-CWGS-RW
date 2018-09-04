#
# python3 SWATH2/SWATH2_COMPUTER.py
#

import os

def replace(x, y, z):
    return x.replace(y, z)

with open("/tmp/A.txt", "r") as w:
    for samp in w:
        samp = samp.strip()
        samp = replace(samp, "-vector", "")
        samp = replace(samp, "-sc1", "")
        samp = replace(samp, "-sc2", "")
        samp = replace(samp, "-sc3", "")
        samp = replace(samp, "-sc4", "")
        samp = replace(samp, "/", "_")
        samp = replace(samp, ".", "_")
        samp = replace(samp, " ", "_")
        samp = replace(samp, "-shATRX-1", "")
        samp = replace(samp, "-shATRX-2", "")
        samp = replace(samp, "-shATRX-3", "")
        samp = replace(samp, "-shATRX-4", "")
        samp = replace(samp, "JFCF-6", "JFCF_6")
        samp = replace(samp, "-shDAXX", "")
        samp = replace(samp, "-vector", "")
        samp = replace(samp, "GM847DM", "GM847")
        samp = replace(samp, "IIIFC", "IIICF")
        samp = replace(samp, "_A", "")
        samp = replace(samp, "_B", "")
        samp = replace(samp, "_C", "")
        samp = replace(samp, "IIICF(P7)", "IIICF_P7")

        print(samp)        