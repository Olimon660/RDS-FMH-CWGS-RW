#
# Apply clustering for VCF files
#
#   python3 cluster.py clusters A.vcf.obj
#

import sys
import vcf
import pickle

# Where we will read our VCF objects
path = sys.argv[1]

# Reference SNPs
ref = pickle.load(sys.argv[2])


