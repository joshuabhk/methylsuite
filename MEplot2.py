"""
Created on Aug 29, 2014

@author: Bong-Hyun Kim
"""
import os
import sys
from sys import stderr

from pysam import Samfile

from seq.genome import Genome2
from call.meplot import GenomeMethylationCounter


usage = '''
%s [-c <chromosome name>|-g <chromosome fasta file>|-o <output_fn>|-s <snp_correction>|-r <trim>|-t|-h] <sam file>

This program will count methylation in the given SAM file.
The SAM file can be output of any generic aligner but specifically tested 
on GSNAP, BSMAP and Bismark.

Options:
-g <chromosome/genome FASTA file>
    FASTA file containing the chromosome reference sequences.
    
-s <snp_correction>    
    Select type of SNP correction. Default value is No correction.
    -s ctsnp for using reverse strand G/(A+G) ratio like in BSMAP
    -s allnt using reverse strand A ratio

-l <trim> (default : 0)
    # of nucleotides to trim from 5' end of reads.

-r <trim> (default : 0)
    # of nucleotide to trim while counting the methylations.
    This option is removing 3' biased region of reads.

-b <quality_base> (default: 33)
    #important for -q and -c option to interprete the quality scores.

-q <quality_trim_cutoff> (default : off)
    This option determines the trim position after appearance of this score.
    
-x <quality_score_exclustion_cutoff> (default : 0)
    Lower than this threshold is not included in methylation counting.

-m <major allele cutoff> (default: 0.5)

-d <major allele detection depth cutoff> (default: 10)

-S
    Include to count reads that are mapped but not the mates are not mapped.
-p
    Print text file for a M-bias plot
-t
    Terse output
-v 
    Debugging outputs via STDERR.
-u
    Including duplicate reads for methylation rates
-h
    Print help message
''' % os.path.basename(sys.argv[0])

if __name__ == '__main__':
    from getopt import getopt

    options, args = getopt(sys.argv[1:], "m:d:r:g:s:o:q:x:b:l:Supvth")
    if not args:
        print usage
        sys.exit()

    samfn = args[0]

    snp_type = 0
    comprehensive_output = True
    output_fn = ''
    output_fp = sys.stdout
    genome_fn = '/scratch/kimb/chr21.fa'
    major_allele_cutoff = 0.5
    major_allele_detection_min_depth = 10
    rtrim = 0
    ltrim = 0
    quality_trim = None
    quality_base = 33
    quality_score_cutoff = None
    mbiasplot_flag = False
    mcount_flag = True
    verbose = 0
    duplicate_including = False
    #all_including = False
    singleton_including = False

    for opt, val in options:
        if opt == '-s':
            if val == 'ctsnp':
                snp_type = 1
            elif val == 'allnt':
                snp_type = 2

        elif opt == '-m':
            major_allele_cutoff = float(val)
        elif opt == '-d':
            major_allele_detection_min_depth = int(val)
        elif opt == '-r':
            rtrim = int(val)
        elif opt == '-l':
            ltrim = int(val)
        elif opt == '-g':
            genome_fn = val
        elif opt == '-t':
            comprehensive_output = False
        elif opt == '-q':
            quality_trim = int(val)
        elif opt == '-b':
            quality_base = int(val)
        elif opt == '-x':
            quality_score_cutoff = int(val)
        elif opt == '-p':
            mbiasplot_flag = True
            mcount_flag = False
	#elif opt == '-A':
        #    all_including = True
	elif opt == '-S':
            singleton_including = True

        elif opt == '-o':
            output_fn = val
            if os.path.exists(output_fn):
                raise Exception("%s already exists" % output_fn)

        elif opt == '-h':
            print usage
            sys.exit()
        elif opt == '-v':
            verbose += 1
        elif opt == '-u':
            duplicate_including = True

    if verbose:
        print >> stderr, "Reading Genome file", genome_fn
    genome = Genome2(genome_fn, verbose=verbose)

    if verbose:
        print >> stderr, "Opening Sam/Bam file", samfn

    samfile = Samfile(samfn)

    methylcounter = GenomeMethylationCounter(genome, samfile,
                                             snp_major_allele_cutoff=major_allele_cutoff,
                                             snp_detection_min_depth=major_allele_detection_min_depth,
                                             rtrim=rtrim,
                                             quality_trim=quality_trim,
                                             quality_base=quality_base,
                                             quality_score_cutoff=quality_score_cutoff,
                                             mcount_flag=mcount_flag,
                                             mbiasplot_flag=mbiasplot_flag,
                                             ltrim=ltrim, verbose=verbose,
                                             duplicate_including=duplicate_including,
                                             singletone_including=singleton_include)
                                             #all_including=all_including)

    if output_fn:
        output_fp = open(output_fn, 'w')

    if mcount_flag:
        methylcounter.build_methylation_bedgraph(sfp=output_fp, snp=snp_type, comprehensive_output=comprehensive_output)

    if mbiasplot_flag:
        methylcounter.build_mbiasplot_output(sfp=output_fp)
