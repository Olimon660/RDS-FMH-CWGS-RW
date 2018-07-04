#
# python STRUCTURAL/STRUCTURAL_SOMATIC.py STAGE2/case_control_input_files_2.txt
#

import os
import sys

def GRIDSS(normal, tumor, out, ass, mod):
    str = 'java -ea -Xmx128g -Dsamjdk.create_index=true -Dgridss.defensiveGC=true -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=true -Dsamjdk.compression_level=1 -cp gridss-1.7.2-gridss-jar-with-dependencies.jar gridss.CallVariants TMP_DIR=8 WORKER_THREADS=1 WORKING_DIR=8 REFERENCE_SEQUENCE=combined.fasta INPUT=%s INPUT=%s OUTPUT=%s ASSEMBLY=%s 2>&1 | tee -a 8/gridss.somatic.%s.log\n'
    return (str % (normal, tumor, out, ass, mod))

with open(sys.argv[1]) as r:
    for line in r:
	    toks = line.strip().split('\t') 
        assert(len(toks) == 2)
        
        #
        # The file follows Erdahl's Excel. IMMORTAL is on first column and is tumor, while MORTAL is normal on the
        # second column.
        #
        
        tumor  = toks[0]
        normal = toks[1]
        
        NORMAL   = '5/REALIGNED_RG_DEDUP_SORTED_HG19_' + normal + '.bam'
        TUMOR    = '5/REALIGNED_RG_DEDUP_SORTED_HG19_' + tumor  + '.bam'
        OUTPUT   = '8/SOMATIC_GRIDSS_' + tumor + '.vcf'
        ASSEMBLY = '8/SOMATIC_GRIDSS_' + tumor + '.bam'
        
        assert(os.path.exists(NORMAL))
        assert(os.path.exists(TUMOR))
        
        cmd = GRIDSS(SAMPLE, OUTPUT, ASSEMBLY, sample)        

        print(cmd)
        os.system(cmd)
