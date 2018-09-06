#
# python3 SWATH2/SWATH2_COMPUTER.py
#

import os

def replace(x, y, z):
    return x.replace(y, z)

with open("/tmp/A.txt", "r") as w:
    for samp in w:
        samp = samp.strip()
        
        if samp == "JFCF-6/T.1/P":
            print("JFCF_6_T_1_P_TEL")
            continue
        elif samp == "JFCF-6/T.1/P (ALT)":
            print("JFCF_6_T_1_P_ALT")
            continue
        elif samp == "IIICF E6E7/C4":
            print("IIICF_E6E7_C4_post")
            continue
        elif samp == "IVG-BF-LXSN_A" or samp == "IVG-BF-LXSN_B" or samp == "IVG-BF-LXSN_C":
            print("IVG_BF_LXSN_pre")
            continue
        elif samp == "LFS05F24 pre_A" or samp == "LFS05F24 pre_B" or samp == "LFS05F24 pre_C":
            print("LFS_05F_24_pre")
            continue
        elif samp == "Met-4A pre_A" or samp == "Met-4A pre_B" or samp == "Met-4A pre_C":
            print("MeT_4A_pre")
            continue
        elif samp == "MeT-4A":
            print("MeT_4A_post")
            continue
        elif samp == "WI38_A" or samp == "WI38_B" or samp == "WI38_C":
            print("WI38")
            continue
        elif samp == "WI-38 VA13/2RA":
            print("VA13")
            continue

        samp = replace(samp, "_A", "")
        samp = replace(samp, "_B", "")
        samp = replace(samp, "_C", "")
        samp = replace(samp, "-", "_")
        samp = replace(samp, "_vector", "")
        samp = replace(samp, "_sc1", "")
        samp = replace(samp, "_sc2", "")
        samp = replace(samp, "_sc3", "")
        samp = replace(samp, "_sc4", "")
        samp = replace(samp, "/", "_")
        samp = replace(samp, ".", "_")
        samp = replace(samp, " ", "_")
        samp = replace(samp, "-shATRX-1", "")
        samp = replace(samp, "-shATRX-2", "")
        samp = replace(samp, "-shATRX-3", "")
        samp = replace(samp, "-shATRX-4", "")
        samp = replace(samp, "JFCF-6", "JFCF_6")
        samp = replace(samp, "_shDAXX", "")
        samp = replace(samp, "_shATRX_1", "")
        samp = replace(samp, "_shATRX_2", "")
        samp = replace(samp, "_shATRX_3", "")
        samp = replace(samp, "_shATRX_4", "")
        samp = replace(samp, "GM847DM", "GM847")
        samp = replace(samp, "IIIFC", "IIICF")
        samp = replace(samp, "IIICF(P7)", "IIICF_P7")
        samp = replace(samp, "IIICF_402_DE_D2", "IIICF_402DE_D2")

        print(samp)        