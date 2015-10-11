#!/usr/bin/env python
'''
Created on Sep 18, 2014

@author: Bong-Hyun Kim
'''
from sys import argv, stderr
from os.path import basename, exists, join
from os import makedirs
from seq.genome import Genome2
from re import sub

verbose = False
skip_errors = True

if __name__ == '__main__':
    usage = '''
    %s <genome_fasta> <output_bed_directory>
    '''%basename(argv[0])
    
    genome_fn = argv[1]
    out_dir = argv[2]
    
    if not exists(out_dir) :
        makedirs(out_dir)
    
    genome = Genome2( genome_fn, verbose=verbose )
    
    for chrname in genome.chromosomes :
        chromosome = genome[chrname]
        
        bed_fn = join(out_dir, chrname + ".methylcontext.bed")
        if exists(bed_fn) :
            if skip_errors:
                print >>stderr, "WARNING: File", bed_fn, "already exists!"
                continue
            else:
                raise Exception( "File aleady exists!", bed_fn )
        
        fp = open(bed_fn, 'w')
        count = 0
        for pos, strand, context in chromosome.itermethylitems() :
            if verbose :
                print >>stderr, chrname, pos, context, strand
            print >>fp, '\t'.join([ chrname, str(pos), str(pos+1), context, strand ])
            count += 1
        fp.close()
        
        if verbose :
            print >>stderr, 'Chromosome', chrname, "has been processed."
            print >>stderr, count, 'records has been written in the bed file', bed_fn + "."
    
    for chrname in genome.chromosomes:
        bed_fn = join(out_dir, chrname + ".methylcontext.bed")
        if not exists(bed_fn):
            print >>stderr, "Error: File", bed_fn, "does not exists!"
            continue
        
        cpg_bed_fn = sub( "methylcontext\.bed", "cpg.bed", bed_fn )
        chg_bed_fn = sub( "methylcontext\.bed", "chg.bed", bed_fn )
        chh_bed_fn = sub( "methylcontext\.bed", "chh.bed", bed_fn )
        
        if exists(cpg_bed_fn):
            print >>stderr, "Error: File", cpg_bed_fn, "already exists!"
            continue
        elif exists(chg_bed_fn):
            print >>stderr, "Error: File", chg_bed_fn, "already exists!"
            continue
        elif exists(chh_bed_fn):
            print >>stderr, "Error: File", chh_bed_fn, "already exists!"
            continue
    
        bed_fp = open( bed_fn )
            
        cpg_fp = open(cpg_bed_fn, 'w')
        chg_fp = open(chg_bed_fn, 'w')
        chh_fp = open(chh_bed_fn, 'w')
        
        for l in bed_fp :
            if "CpG" in l:
                cpg_fp.write(l)
            elif "CHG" in l:
                chg_fp.write(l)
            elif "CHH" in l:
                chh_fp.write(l)
            else :
                print >>stderr, "WARNING: Unexpected Line format!"
                stderr.write(l)
                
        cpg_fp.close()
        chg_fp.close()
        chh_fp.close()
        