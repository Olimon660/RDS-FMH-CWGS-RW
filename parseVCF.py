#
# python3 parseVCF.py vcfs A.vcf 6
#

import os
import sys
import vcf
import pickle
from os import listdir
from os.path import isfile, join

def getKey(x):
    return str(x.CHROM) + ':' + str(x.POS)

def readVCF(file, ref = None):
    print('Reading: ' + file)
    x = {}
    for r in vcf.Reader(open(file)):
        key = getKey(r)        
        if ref is None:
            x[key] = { 'CHROM':r.CHROM, 'POS':r.POS, 'REF':r.REF, 'ALT':r.ALT, 'key':key }
        else:
            assert(len(r.samples) == 1)            
            if key in ref:            
                x[key] = { 'CHROM':r.CHROM, 'POS':r.POS, 'REF':r.REF, 'ALT':r.ALT, 'key':key, 'GT':r.samples[0].gt_type }
    return x

def writeVCF(path, file, x):
    with open(path + os.sep + os.path.basename(file) + '.obj', 'wb') as w:
        pickle.dump(list(x), w, pickle.HIGHEST_PROTOCOL)    

# Read common reference SNPs
ref = readVCF(sys.argv[2])

# Write reference SNPs
writeVCF(sys.argv[1], sys.argv[2], ref)

# Read all samples
samps = [f for f in listdir(sys.argv[3]) if isfile(join(sys.argv[3], f))]

for samp in samps:
    if samp.endswith('.vcf'):
        vcf = readVCF(sys.argv[3] + os.sep + samp, ref)
        print(len(vcf))
        writeVCF(sys.argv[1], sys.argv[3] + os.sep + samp, vcf)
