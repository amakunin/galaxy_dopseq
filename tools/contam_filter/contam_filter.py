#!/usr/bin/env python

import sys
import tempfile
import os
import pysam


# minimum MAPQ of the read retained after contamination filtering
min_qual = int(sys.argv[1])

# inputs: bams vs target and contam genomes
t_filename = sys.argv[2]
c_filename = sys.argv[3]

# output: filtered bam alignment
filter_filename = sys.argv[4]
#contam_filename = sys.argv[5]
#unmap_filename = sys.argv[6]

# sort inputs by read name in temp
t = tempfile.NamedTemporaryFile(suffix = '_t.bam', delete=False)
c = tempfile.NamedTemporaryFile(suffix = '_c.bam', delete=False)
tname = t.name
cname = c.name
t.close()
c.close()
#pysam.sort('-n', '-o', tname, t_filename) # samtools 0.1.19
#pysam.sort('-n', '-o', cname, c_filename) # samtools 0.1.19
pysam.sort('-n', t_filename, tname[:-4])
pysam.sort('-n', c_filename, cname[:-4])

# create intermediate filtered output
f = tempfile.NamedTemporaryFile(suffix = '_f.bam', delete=False)
fname = f.name
f.close()

with pysam.AlignmentFile(tname) as tfile, pysam.AlignmentFile(cname) as cfile:
    filter_file = pysam.AlignmentFile(fname,'wb', template=tfile)
    #contam_file = pysam.AlignmentFile(contam_filename,'wb', template=cfile)
    #unmap_file = pysam.AlignmentFile(unmap_filename,'wb', template=tfile)
    i = 0
    for tread in tfile:
        cread = cfile.next()
        assert tread.query_name == cread.query_name
        if tread.mapping_quality >= min_qual and tread.mapping_quality >= cread.mapping_quality:
            i += 1
            filter_file.write(tread)
        #if tread.mapping_quality < min_qual: # unmapped
        #    unmap_file.write(tread) # output mapping for target genome
        #elif tread.mapping_quality < cread.mapping_quality: # contamination
        #    contam_file.write(cread) # output mapping for contamination genome
        #else: #target
        #    filter_file.write(tread)
    filter_file.close()
    #contam_file.close()
    #unmap_file.close()

# sort 
pysam.sort(fname, fname[:-4]+'_srt')
os.rename(fname[:-4]+'_srt.bam', filter_filename)
sys.stdout.write('Writing %d alignments to filtered output file %s\n'%(i, filter_filename))
#sys.stdout.write('3 output files: filtered %s, contamination %s, unmapped %s\n'%(filter_filename, contam_filename, unmap_filename))

#clean-up
os.unlink(tname)
os.unlink(cname)
os.unlink(fname)