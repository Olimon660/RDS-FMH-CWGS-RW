#
# python3 SWATH/SWATH_CLEANUP.py
#

import os

frds = [ "IIICF", "IIICF_d2", "IIICF_402DE_D2", "IIICF_E6_A1", "IIICF_E6_A2", "IIICF_E6E7_A1", "IIICF_T_B1", "IIICF_T_B3",   \
         "IIICF_T_C3", "IIICF_T_C4", "JFCF_6_T_1_C", "JFCF_6_T_1_D", "JFCF_6_T_1_F", "JFCF_6_T_1_G", "JFCF_6_T_1_H",         \
         "JFCF_6_T_1_L", "JFCF_6_T_1_P_TEL", "JFCF_6_T_1_P_ALT", "JFCF_6_T_1_R", "JFCF_6_T_1J_11C", "JFCF_6_T_1J_11E",       \
         "JFCF_6_T_2H", "JFCF_6_T_5K", "GM02063", "GM847", "IIICF_E6E7_C4_pre", "IIICF_E6E7_C4_post", "IIICF_a2", "IIICF_c", \
         "IIICF_T_A6", "IVG_BF_LXSN_pre", "IVG_BF_LXSN_post", "JFCF_6", "JFCF_6_T_1_M", "JFCF_6_P_pLKO_5", "JFCF_6_T_1_Q",   \
         "JFCF_6_T_1J_1_3C", "JFCF_6_T_1J_6B", "LFS_05F_24_post", "LFS_05F_24_pre", "MeT_4A_post", "MeT_4A_pre", "WI38",  \
         "VA13", "IIICF_P7", "IIICF_P9" ]

def noRep(x):
    return x.replace("_r1", "").replace("_r2", "").replace("_r3", "").replace("_r4", "")

def onlyRep(x):
    return x.split("_")[-1]

def findFri(x):
    assert("_r1" in x or "_r2" in x or "_r3" in x or "_r4" in x)
    rp = onlyRep(x)

    x = noRep(x)    
    if x == "MeT-4A":
        x = "MeT_4A_post"        
    x = x.replace("/", "_")
    x = x.replace(" ", "_")    
    x = x.replace("GM847DM", "GM847")
    x = x.replace("IIICF-T", "IIICF_T")
    x = x.replace("(P7)", "_P7")
    x = x.replace("(P9)", "_P9")
    x = x.replace("IIIFC", "IIICF")
    x = x.replace("IVG-BF", "IVG_BF")
    x = x.replace("-LXSN", "_LXSN")
    x = x.replace("JFCF-6", "JFCF_6")
    x = x.replace(".1", "_1")
    x = x.replace(".2", "_2")
    x = x.replace(".3", "_3")
    x = x.replace(".4", "_4")
    x = x.replace(".5", "_5")
    x = x.replace("(ALT)", "ALT")
    x = x.replace("-sc2", "_sc2")
    x = x.replace("-", "_")
    x = x.replace("LFS05F24", "LFS_05F_24")
    x = x.replace("Met", "MeT")

    if x == "JFCF_6_T_1_P":
       x = "JFCF_6_T_1_P_TEL"
    elif x == "IIICF_E6E7_C4":
        x = "IIICF_E6E7_A1"
    elif x == "IVG_BF_LXSN":
        x = "IVG_BF_LXSN_pre"
    elif x == "JFCF_6_T_1_P_sc2":
        x = "JFCF_6_P_pLKO_5"
    elif x == "WI_38_VA13_2RA":
        x = "VA13"

    assert(x in frds)
    return x + "_" + rp

cns = {}
new = []

with open("SWATH/sample.txt", "r") as r:
    for line in r:
        samp = line.strip()
        
        if not samp in cns:
            cns[samp] = 1
        else:
            cns[samp] = cns[samp] + 1
        
        if samp == "IIICF/a2_A":
            new.append("IIICF/a2_r4")            
        elif samp.endswith("_A"):
            new.append(samp.replace("_A", "_r1"))
        elif samp.endswith("_B"):
            new.append(samp.replace("_B", "_r2"))
        elif samp.endswith("_C"):
            new.append(samp.replace("_C", "_r3"))
        else:
            new.append(samp + "_r" + str(cns[samp]))

with open("SWATH/newSample.txt", "w") as w:
    for samp in new:
        # Make sure it's a user friendly name
        w.write(findFri(samp) + "\n")
        