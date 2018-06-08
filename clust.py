import vcf

def readVCF(file):
    for rec in vcf.Reader(open(file)):
        print(rec)

# Read common SNPS        
readVCF('CEU.exon.2010_03.genotypes.vcf')

# Read all samples

# Genereate a distance matrix
