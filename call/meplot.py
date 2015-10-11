'''
Created on Aug 27, 2014

@author: Bong-Hyun Kim
'''

from array import array
from cStringIO import StringIO
from seq.fastq import quality_string2score
from tempfile import NamedTemporaryFile
import h5py
import os, time
from sys import stderr, exit
from util import overlap, split_middle, RangeInclusionError

# import sys
# arraytype = 'B' #1 byte range: 0 - 255
# arraytype = 'I'  # 2 bytes range: 0 - 2^16-1

class NotMatchingChromosomeError( Exception ):
    pass


class InvalidChromosomeNameSamfileError(Exception):
    pass


class GenomeMethylationCounter :
    def __init__(self, genome, samfile, **kwargs):
        
        # Assume that the samfile.references 
        # are in the same as the sorted chromosome order
        # In fact the samfile should be ordered by coordinates
        # or at least sorted by chromosome.
        self.samfile = samfile
        self.kwargs = kwargs
        self.counters = []
        
        #assign kwargs to self 
        for k,v in kwargs.iteritems() :
            setattr(self, k, v)
        
        if self.verbose :
            print >>stderr, time.ctime()
            print >>stderr, "Arguments passed from Main function:"
            print >>stderr, kwargs
        
        for chromosome_name in genome.chromosome_names :
            if self.verbose :
                print >>stderr, "Setting up Chromosome", chromosome_name
            chromosome = genome[ chromosome_name ]
            
            if self.verbose :                    
                print >>stderr, "Counting Chromosome", chromosome_name
            counter = ChromosomeMethylationCounter( chromosome, samfile, **kwargs )
            self.counters.append( counter )
            
            if self.verbose:
                print >>stderr, "Done", chromosome_name, time.ctime() 
            continue
            
    
    def build_mbiasplot_output(self, sfp=None):
        for counter in self.counters :
            counter.build_mbiasplot_output(sfp=sfp)
    
    def build_methylation_bedgraph(self, sfp=None, snp=0, comprehensive_output=True):
        for counter in self.counters :
            counter.build_methylation_bedgraph( sfp=sfp, 
                                                snp=snp, 
                                                comprehensive_output=comprehensive_output)
    def build_methylation_summary(self, sfp=None ):
        for counter in self.counters:
            counter.build_methylation_summary(sfp=sfp)
                
class AmbiguousReadPairError(Exception):
    pass


class Read1Read2SettingError(Exception):
    pass

class ChromosomeIndexOutOfBoundError(Exception):
    pass


class RawReadIndexOutOfBoundError(Exception):
    pass


class ChromosomeMethylationCounter:

    def __init__(self, chromosome, samfile,
                 snp_major_allele_cutoff=0.5, snp_detection_min_depth=10,
                 rtrim=0, quality_trim=None, quality_base=33, 
                 quality_score_cutoff=None,
                 autostart=True, readlength=200, mcount_flag=True,
                 mbiasplot_flag=False, ltrim=0, 
                 arraytype='I', tempdir='./', tempfilename='', suffix='.hdf5',
                 count_only_paired_reads=True, verbose=False,
                 duplicate_including=False,
                 #all_including=False,
                 singleton_including=False,
                 trim_type='simple'
                 ):
        
        self.chromosome = chromosome
        self.name = chromosome.name
        self.get_sequence = chromosome.get_sequence  # for faster access
        self.get_nucleotide = chromosome.get_nucleotide  # for faster access
        self.get_methylation_context = chromosome.get_methylation_context
        self.iterposnt = chromosome.iterposnt
        self.itermethylposseq = chromosome.itermethylposseq
        self.itermethylitems = chromosome.itermethylitems
        
        self.snp_major_allele_cutoff = snp_major_allele_cutoff
        self.snp_detection_min_depth = snp_detection_min_depth
        self.quality_trim = quality_trim
        self.quality_base = quality_base
        self.quality_score_cutoff = quality_score_cutoff

        self.readlength = readlength
        self.mcount_flag = mcount_flag
        self.mbiasplot_flag = mbiasplot_flag

        self.rtrim = rtrim
        self.ltrim = ltrim
        self.arraytype=arraytype
        self.tempdir = tempdir
        self.tempfilename=tempfilename
        self.suffix=suffix
        self.count_only_paired_reads = count_only_paired_reads
        self.verbose=verbose
        self.duplicate_including=duplicate_including
        self.singleton_including=singleton_including
        #self.all_including=all_including
        self.delayed_analysis_pool = {}
        
        self.progress = -1 #progress indicator for verbose output
        self.trim_type = trim_type
        
        if not self.tempfilename :
            self.tempfile_fp = NamedTemporaryFile(suffix=suffix, dir=tempdir)
            self.tempfilename = self.tempfile_fp.name
        else :
            if os.path.exists( tempfilename ) :
                raise Exception( "File already exists!", tempfilename )

        if self.mcount_flag:
            self.length = chromosome.get_length()
            self.counts = self.initialize_counts()
            
        if self.mbiasplot_flag:
            self.readcounts_contexts = self.initialize_mplot_counts()

        self.samfile = samfile
        if samfile and autostart and (self.mcount_flag or self.mbiasplot_flag):
            analyze = self.analyze_read
            self.tid = samfile.gettid(self.name)

            if self.tid == -1 and self.mcount_flag == True:
                #print >>sys.stderr, "Chromosome name is not valid!"
                #sys.exit(-1)
                raise InvalidChromosomeNameSamfileError(self.samfile, self.tid)

            #for aread in samfile:
                #analyze(aread)
            #the samfile needs to be indexed.
            for aread in samfile.fetch( reference=self.chromosome.name ) : #, callback=analyze )
                analyze(aread)
                        
        #hdf5 saving routine.
        #not that useful!!
        #if self.mcount_flag :
            #self.counts = self.create_hdf5_datasets()

    def initialize_mplot_counts(self):
        readlength = self.readlength
        readcounts_contexts = {'CpG':{
                '+':{'M':array('L', [0]) * readlength, 
                    'U':array('L', [0]) * readlength, 
                    'O':array('L', [0]) * readlength}, 
                '-':{'M':array('L', [0]) * readlength, 
                    'U':array('L', [0]) * readlength, 
                    'O':array('L', [0]) * readlength}}, 
            'CHG':{
                '+':{'M':array('L', [0]) * readlength, 
                    'U':array('L', [0]) * readlength, 
                    'O':array('L', [0]) * readlength}, 
                '-':{'M':array('L', [0]) * readlength, 
                    'U':array('L', [0]) * readlength, 
                    'O':array('L', [0]) * readlength}}, 
            'CHH':{
                '+':{'M':array('L', [0]) * readlength, 
                    'U':array('L', [0]) * readlength, 
                    'O':array('L', [0]) * readlength}, 
                '-':{'M':array('L', [0]) * readlength, 
                    'U':array('L', [0]) * readlength, 
                    'O':array('L', [0]) * readlength}}}
        return readcounts_contexts

    def initialize_counts(self):
        length = self.length
        arraytype = self.arraytype
        counts = {'+':{'C':array(arraytype, [0]) * length, 'G':array(arraytype, [0]) * length, 
                'A':array(arraytype, [0]) * length, 
                'T':array(arraytype, [0]) * length,
                'N':array(arraytype,[0]) * length }, 
            '-':{'C':array(arraytype, [0]) * length, 
                'G':array(arraytype, [0]) * length, 
                'A':array(arraytype, [0]) * length, 
                'T':array(arraytype, [0]) * length,
                'N':array(arraytype, [0]) * length }
            }
        return counts

    def create_hdf5_datasets(self):
        tempfilename = self.tempfilename
        arraytype = self.arraytype
        length = self.length
        
        counts = h5py.File(tempfilename, mode='w')
        counts.create_dataset("+/C", (length, ), arraytype, data=self.counts['+']['C'])
        counts.create_dataset("+/G", (length, ), arraytype, data=self.counts['+']['G'])
        counts.create_dataset("+/A", (length, ), arraytype, data=self.counts['+']['A'])
        counts.create_dataset("+/T", (length, ), arraytype, data=self.counts['+']['T'])
        counts.create_dataset("+/N", (length, ), arraytype, data=self.counts['+']['N'])
        counts.create_dataset("-/C", (length, ), arraytype, data=self.counts['-']['C'])
        counts.create_dataset("-/G", (length, ), arraytype, data=self.counts['-']['G'])
        counts.create_dataset("-/A", (length, ), arraytype, data=self.counts['-']['A'])
        counts.create_dataset("-/T", (length, ), arraytype, data=self.counts['-']['T'])
        counts.create_dataset("-/N", (length, ), arraytype, data=self.counts['-']['N'])
        
        #clear the contents.
        del self.counts['+']['C'][:]
        del self.counts['+']['G'][:]
        del self.counts['+']['A'][:]
        del self.counts['+']['T'][:]
        del self.counts['+']['N'][:]
        
        del self.counts['-']['C'][:]
        del self.counts['-']['G'][:]
        del self.counts['-']['A'][:]
        del self.counts['-']['T'][:]
        del self.counts['-']['N'][:]
        
        self.counts['-'].clear()
        self.counts['+'].clear()
        self.counts.clear()
        
        return counts

        
    def analyze_read(self, aread):
        # callback for the Samfile fetch function
        # nothing to do if the read not mapped
        if aread.is_unmapped:
            return

        #properly paired read will be read 
        #if count_only_paired_reads option is true
        if aread.is_paired and self.count_only_paired_reads :
            if not aread.is_proper_pair :
		if aread.mate_is_unmapped and self.singleton_including :
			pass
		else :
                	return            
        
        if aread.is_duplicate and not self.duplicate_including:
            return
        
        # when aligned to the correct reference chromosome
        if self.tid == -1 or self.tid == aread.tid:
            # skip multimappers
            if aread.is_secondary:
                return
            try:
                if aread.opt('NH') > 1:
                    # print aread
                    return
            except KeyError:
                pass

            '''
            refseq, queryseq, qualseq = self.get_aligned_sequences( aread )
            if '-' in refseq or '-' in queryseq :
                print aread.aligned_pairs
                print aread.qstart, aread.qend
                print refseq
                print queryseq
                print qualseq
                refseq, queryseq, qualseq = self.get_aligned_sequences_simple( aread )
                print refseq
                print queryseq
                print qualseq
                print "\n"
            refseq2, queryseq2, qualseq2 = self.get_aligned_sequences_simple( aread )
            assert( refseq == refseq2 )
            assert( queryseq == queryseq2 )
            assert( qualseq == qualseq2 )
            '''

            if not aread.is_paired or not aread.is_proper_pair : #single or non-proper pair cases
                refseq, queryseq, qualseq = self.get_aligned_sequences_simple(aread)
                qualscore = quality_string2score(qualseq)
    
                if self.mcount_flag:
                    self.count_methylation(aread, refseq, queryseq, qualscore)
    
                if self.mbiasplot_flag:
                    self.count_mbiasplot(aread, refseq, queryseq, qualscore)
                
                if self.verbose :
                    progress = int(aread.aend*100.0/self.length)
                    if progress > self.progress :
                        self.progress = progress
                        print >>stderr, time.ctime(), self.chromosome.name, '%d%%'%self.progress
            else : #analyze properly paired reads
                if aread.pnext > aread.pos :
                    self.add_delayed_analysis_pool(aread)
                elif aread.pnext < aread.pos :
                    aread1 = self.get_left_read( aread )
                    
                    if not aread1 :
                        print >>stderr, "Error: aread1 not found!", aread.qname
                        if aread.is_proper_pair : print >>stderr, 'properly paired'
                        print >>stderr, 'pos:', aread.pos, 'pnext:', aread.pnext, 
                        print >>stderr, 'paired:', aread.is_paired, 'read1:', aread.is_read1
                        print >>stderr, 'seq:', aread.seq
                        
                        return
                    
                    #test if the read1 and read2 is set correctly
                    if aread1.is_read1 and aread.is_read2 :
                        pass
                    elif aread1.is_read2 and aread.is_read1 :
                        #need to swap!
                        aread, aread1 = aread1, aread
                    else :
                        #raise Read1Read2SettingError( aread, aread1 )
			print >>stderr, "Error: aread1 and not aread2 are not matching", aread.qname, aread1.qname
			print >>stderr, aread
                        if aread.is_proper_pair : print >>stderr, 'properly paired'
                        print >>stderr, 'seq:', aread.seq
			print >>stderr, 'qual:', aread.qual
			print >>stderr, aread1
                        if aread1.is_proper_pair : print >>stderr, 'properly paired'
                        print >>stderr, 'seq:', aread1.seq
			print >>stderr, 'qual:', aread1.qual
                        
                        return
                    

                    
                    self.analyze_paired_read(aread1, aread)
                    
                elif aread.pnext == aread.pos :
                    aread1 = self.get_left_read( aread )
                    if aread1 :
                        self.analyze_paired_read(aread1, aread)
                    else :
                        self.add_delayed_analysis_pool(aread)
                    
        else :
            #When the alignment is for the different reference genome
            raise NotMatchingChromosomeError( self.name )

    def add_delayed_analysis_pool(self, aread):
        if aread.pos in self.delayed_analysis_pool :
            self.delayed_analysis_pool[aread.pos].append( aread )
        else :
            self.delayed_analysis_pool[aread.pos] = [ aread ]
    
    def get_left_read(self, aread):
        '''
        Find a paired read of the given aread.
        returns AlignedRead object instance if the pair found,
        otherwise returns None.
        '''
        if aread.pnext not in self.delayed_analysis_pool :
            return None
        
        read_candidates = []
        current_pos = aread.pos
        for read in self.delayed_analysis_pool[aread.pnext] :
            if read.pnext == current_pos :
                if read.qname and read.qname == aread.qname :
                    self.delayed_analysis_pool[aread.pnext].remove( read )
                    
                    if not self.delayed_analysis_pool[aread.pnext] :
                        del self.delayed_analysis_pool[aread.pnext]
                    
                    return read
                else :
                    read_candidates.append(read)
                    
        if len(read_candidates) > 1 :
            raise AmbiguousReadPairError( aread, read_candidates )
        elif len(read_candidates) == 1 :
            self.delayed_analysis_pool[aread.pnext].remove(read_candidates[0])
            if not self.delayed_analysis_pool[aread.pnext] :
                del self.delayed_analysis_pool[aread.pnext]
            return read_candidates[0]
        else :
            return None


    def analyze_paired_read(self, aread1, aread2):
        #determine center point between the two reads
        #after post-mortem trimming
        
        if self.mcount_flag:
            iqual1 = quality_string2score(aread1.qual, self.quality_base)
            iqual2 = quality_string2score(aread2.qual, self.quality_base)
            self.count_methylation_directional_pair(aread1, aread2, iqual1, iqual2)
            
        if self.mbiasplot_flag:
            refseq1, queryseq1, qualseq1 = self.get_aligned_sequences_simple(aread1)
            qualscore1 = quality_string2score(qualseq1)
            refseq2, queryseq2, qualseq2 = self.get_aligned_sequences_simple(aread2)
            qualscore2 = quality_string2score(qualseq2)

            self.count_mbiasplot(aread1, refseq1, queryseq1, qualscore1)
            self.count_mbiasplot(aread2, refseq2, queryseq2, qualscore2)
        
        if self.verbose :
            progress = int(aread2.aend*100.0/self.length)
            if progress > self.progress :
                self.progress = progress
                print >>stderr, time.ctime(), self.chromosome.name, '%d%%'%self.progress
        
    def infer_reference_sequence_length_from_cigar(self, cigar):
        length = 0
        for op, le in cigar:
            if op == 0 or op == 2 or op == 4 or op == 7 or op == 8:
                length += le

        return length

    def get_aligned_sequences_simple(self, aread):
        # returns aligned reference sequence
        # this function should be more stable than the get_aligned_sequences 
        # function, since this one uses pre-calculated pysam information effciently,
        # whereas the get_aligned_sequences function uses CIGAR string to infer
        # alignments.
        
        #Question: Is this really necessary?

        query = aread.query
        qqual = aread.qqual

        queseq = []
        quaseq = []
        refseq = []

        #chromosome_name = self.samfile.getrname(aread.tid)

        for qpos, rpos in aread.aligned_pairs:
            if qpos == None:
                queseq.append('-')  # for debugging purpose
                quaseq.append('-')  # for debugging purpose
            else:
                queseq.append(query[qpos])
                quaseq.append(qqual[qpos])

            if rpos == None:
                refseq.append('-')
            else:
                refseq.append(self.get_nucleotide(rpos)) #chromosome_name, rpos))

        #if self.verbose :
            #print >>stderr, "refseq:", refseq
            #print >>stderr, 'queseq:', queseq
            #print >>stderr, 'quaseq:', quaseq
            
        return ''.join(refseq).upper(), ''.join(queseq).upper(), ''.join(quaseq)

    def get_aligned_sequences(self, aread):
        # returns aligned reference sequence
        # using starting position and CIGAR operation
        # Notably, the sequence is complemented compared to the + strand
        #(not reverse complement)
        # if the read is mapper to the reverse strand (or - strand).

        qseq = aread.query
        qual = aread.qqual
        refstart = aread.pos

        # reflen seems to be excluding soft clipped region
        reflen = aread.alen
        #reflen = self.infer_reference_sequence_length_from_cigar( aread.cigar )

        # SAM format always represent on the forward strand!
        strand = '+'
        #chromosome_name = self.samfile.getrname(aread.tid)
        # returns both sequence and direction
        
        #old function call from Genome class
        #rseq, direction = self.get_sequence(  # @UnusedVariable
        #    chromosome_name, refstart, reflen, strand)
        #new function call from Chromosome class
        rseq, direction = self.get_sequence(refstart, reflen, strand)  # @UnusedVariable

        if not rseq.isupper():
            rseq = rseq.upper()
        if not qseq.isupper():
            qseq = qseq.upper()

        rpos = 0
        qpos = 0
        refseq = []
        queseq = []
        qualseq = []
        for operation, length in aread.cigar:
            # match 'M' or seqequal '='
            if operation == 0 or operation == 7 or operation == 8:
                refseq.append(rseq[rpos:rpos + length])
                queseq.append(qseq[qpos:qpos + length])
                qualseq.append(qual[qpos:qpos + length])
                rpos += length
                qpos += length
            # insertion 'I'
            elif operation == 1:
                refseq.append('-' * length)
                queseq.append(qseq[qpos:qpos + length])
                qualseq.append(qual[qpos:qpos + length])
                qpos += length
            # deletion D or reference skip N
            elif operation == 2 or operation == 3:
                refseq.append(rseq[rpos:rpos + length])
                queseq.append('-' * length)
                qualseq.append('-' * length)
                rpos += length
            # soft clipping S
            elif operation == 4:
                # original idea need to include them
                #refseq += rseq[rpos:rpos+length].lower()
                #queseq += qseq[qpos:qpos+length].lower()
                #rpos += length
                #qpos += length
                # But it might be better to exclude it
                # from the alignemnt due to mapping start and end positions are
                # all without the soft clip
                pass
            # hard clipping H
            elif operation == 5:
                raise Exception("Hard clipping is not defined")
            # padding
            elif operation == 6:
                raise Exception("Padding is not defined")

        try:
            assert(len(refseq) == len(queseq))
            assert(qpos == aread.qlen)
            assert(reflen == rpos)
        except:
            print aread.cigar
            print rseq
            print qseq
            print aread.qstart, aread.qend, aread.qlen
            print 'refseq', refseq
            print 'queseq', queseq
            raise

        # debugging
        '''
        diffcount = sum( [1 for r, q in zip( refseq, queseq ) if r != q] )
        if diffcount*1.0/len(refseq) > 0.7 :
            print >>sys.stderr, "WARNING: reference sequence and query sequence is quite different. Diffcount %:", diffcount*1.0/len(refseq)
            print >>sys.stderr, aread.cigar
            print >>sys.stderr, rseq
            print >>sys.stderr, qseq
            print >>sys.stderr, aread.qstart, aread.qend, aread.qlen
            print >>sys.stderr, 'refseq', refseq
            print >>sys.stderr, 'queseq', queseq
        '''

        # debugging
        '''
        if aread.cigar[0][0] == 4 :
            print >>sys.stderr, aread
            print >>sys.stderr, aread.cigar
            print >>sys.stderr, rseq
            print >>sys.stderr, qseq
            print >>sys.stderr, aread.qstart, aread.qend, aread.qlen
            print >>sys.stderr, 'refseq', refseq
            print >>sys.stderr, 'queseq', queseq
        '''

        # debugging
        '''
        print aread
        print rseq
        print qseq
        print aread.qstart, aread.qend, aread.qlen
        print 'refseq', refseq
        print 'queseq', queseq
        '''

        return ''.join(refseq), ''.join(queseq), ''.join(qualseq)

    def determine_trim(self, aread, quality_score, queryseq):
        # find the trim point according to the class variable
        # quality_trim
        # this function returns the number of nucleotides to be trimmed.

        # print "is_reverse:", aread.is_reverse
        # print "qqual:", aread.qqual
        # print " qual:", aread.qual

        threshold = self.quality_trim

        if aread.is_reverse:
            quality_score = quality_score[::-1]
            queryseq = queryseq[::-1]  # reversing the string!

        trim = 0
        for i, (q, a) in enumerate(zip(quality_score, queryseq)):
            if a.isalpha() and q <= threshold:
                trim = len(quality_score) - i
                break

        # print "qual:", aread.qual
        # print "trim:", trim

        return trim

    def determine_trim_index(self, aread, quality_score, queryseq):
        # find the trim point according to the class variable
        # quality_trim
        # this function returns index of the trim position.
        # The return index is always 5' to 3' in the read,
        # not as in aligned reads in SAM format.

        # print "is_reverse:", aread.is_reverse
        # print "qqual:", aread.qqual
        # print " qual:", aread.qual

        threshold = self.quality_trim

        # if aread.is_reverse :
        #quality_score = quality_score[::-1]
        # queryseq = queryseq[::-1] #reversing the string!

        if not aread.is_reverse:
            #(+) strand mapped read
            offset = aread.qstart
            for i, (q, a) in enumerate(zip(quality_score, queryseq)):
                if a.isalpha() and q <= threshold:
                    return i + offset
            return aread.qend
        else:
            #(-) strand mapped read
            offset = aread.rlen - aread.qend
            for i, (q, a) in enumerate(reversed(zip(quality_score, queryseq))):
                if a.isalpha() and q <= threshold:
                    return i + offset
            return aread.rlen - aread.qstart
        

    def add_integer_quality_scores(self, aread):
        '''
        add integer list of quality scores as iqqual (from qqual)
        and iqual (from qual).
        '''
        if self.verbose :
            print >>stderr, "aread.iqual started to set."
            
        #test if aread has interger transformed qqual
        #if not hasattr(aread, "iqqual"):
        #    aread.iqqual = quality_string2score(aread.qqual, self.quality_base)
        
        #try :
            #if len(aread.qqual) == len(aread.qual):
            #    aread.iqual = aread.iqqual
            #else:
            #    aread.iqual = quality_string2score(aread.qual, self.quality_base)
        #    aread.iqual
        #except AttributeError :
        #    aread.iqual = quality_string2score(aread.qual, self.quality_base)
        
        aread.iqual = quality_string2score(aread.qual, self.quality_base)    
        if self.verbose :
            print >>stderr, "aread.iqual has been set."

    def determine_trimmed_range(self, aread, iqual):
        #is supposedly better and faster than determine_trim_index.
        #returns ranges of query or qqual according to the trim type and the self.trim 
        #the ranges are in RAW READS INDEXED.
                
        trim_type = self.trim_type
        threshold = self.quality_trim
        
        #check the assumption
        assert( aread.qlen == aread.qend - aread.qstart )
        assert( trim_type == 'simple' )# currently only simple trimming is supported!
        
        #self.add_integer_quality_scores(aread)
        iiqual = list(enumerate(iqual)) #indexed iqqual
        
        if aread.is_reverse :
            iiqual.reverse()
        
        #trimming algorithm
        for i, q in iiqual :
            if q <= threshold :
                break
                    
        if not aread.is_reverse :
            start = iiqual[0][0]
            end = i
        else :
            end = iiqual[0][0] + 1
            start = i+1
        
        if self.verbose > 3 :
            print >>stderr, "****"
            print >>stderr, 'determine_trimmedD_range'
            print >>stderr, aread
            print >>stderr, 'iiqual:', iiqual
            print >>stderr, 'start:',start,'end:', end
            print >>stderr, 'threshold:', threshold
            print >>stderr, ""
            print >>stderr, 'return:', max(start, aread.qstart), min(end, aread.qend)
            print >>stderr, '****'
             
        return max(start, aread.qstart), min(end, aread.qend)
    

    def count_mbiasplot(self, aread, refseq, queryseq, quality_score):
        # update the class variables using the aligned refseq and queryseq
        # only counts with upper characters!!

        rtrim = self.rtrim
        ltrim = self.ltrim
        cutoff = self.quality_score_cutoff

        # new way of incorporating ltrim and rtrim (5' and 3' end trimming)
        counting_start = ltrim
        counting_end = aread.rlen - rtrim

        if self.quality_trim != None:
            counting_end = min(
                counting_end, self.determine_trim_index(aread, quality_score, queryseq))

        # forward read
        if not aread.is_reverse:
            cpgm, cpgu, cpgo = self.readcounts_contexts['CpG'][
                '+']['M'], self.readcounts_contexts['CpG']['+']['U'], self.readcounts_contexts['CpG']['+']['O']
            chgm, chgu, chgo = self.readcounts_contexts['CHG'][
                '+']['M'], self.readcounts_contexts['CHG']['+']['U'], self.readcounts_contexts['CHG']['+']['O']
            chhm, chhu, chho = self.readcounts_contexts['CHH'][
                '+']['M'], self.readcounts_contexts['CHH']['+']['U'], self.readcounts_contexts['CHH']['+']['O']
            offset = aread.qstart
            roffset = aread.pos

            # not right!
            # tried to use ltrim but it makes things difficult due to the start position counting.
            # look at the above and blew for codes implementing ltrim and rtrim
            # if rtrim :
            # cut few nt at the 3' terminal
            #refseq = refseq[ltrim:-rtrim]
            #queryseq = queryseq[ltrim:-rtrim]
            #quality_score = quality_score[ltrim:-rtrim]
            # elif ltrim :
            #refseq = refseq[ltrim:]
            #queryseq = queryseq[ltrim:]
            #quality_score = quality_score[ltrim:]

            #chromosome_name = self.samfile.getrname(aread.tid)
            for r, q, s in zip(refseq, queryseq, quality_score):
                if offset >= counting_start:
                    if offset >= counting_end:
                        break
                    else:
                        if q.isalpha() and s >= cutoff:
                            if r == 'C':
                                methylcontext = self.get_methylation_context(
                                                        roffset, direction=1)
                                if methylcontext == 'CHH':
                                    if q == 'C':
                                        chhm[offset] += 1
                                    elif q == 'T':
                                        chhu[offset] += 1
                                    else:
                                        chho[offset] += 1
                                elif methylcontext == 'CHG':
                                    if q == 'C':
                                        chgm[offset] += 1
                                    elif q == 'T':
                                        chgu[offset] += 1
                                    else:
                                        chgo[offset] += 1
                                else:  # methylcontext == 'CpG' :
                                    if q == 'C':
                                        cpgm[offset] += 1
                                    elif q == 'T':
                                        cpgu[offset] += 1
                                    else:
                                        cpgo[offset] += 1

                if q.isalpha():
                    offset += 1
                if r.isalpha():
                    roffset += 1

        # backward strand query
        else:
            cpgm, cpgu, cpgo = self.readcounts_contexts['CpG'][
                '-']['M'], self.readcounts_contexts['CpG']['-']['U'], self.readcounts_contexts['CpG']['-']['O']
            chgm, chgu, chgo = self.readcounts_contexts['CHG'][
                '-']['M'], self.readcounts_contexts['CHG']['-']['U'], self.readcounts_contexts['CHG']['-']['O']
            chhm, chhu, chho = self.readcounts_contexts['CHH'][
                '-']['M'], self.readcounts_contexts['CHH']['-']['U'], self.readcounts_contexts['CHH']['-']['O']
            # read start position for the backward case!
            # Note that the read is reversed!
            offset = aread.rlen - aread.qend  # backward adjustment
            roffset = aread.aend - 1  # adjustment for reversing

            assert(offset >= 0)

            # if ltrim :
            # cut few nt at the 3' terminal
            #refseq = refseq[rtrim:-ltrim]
            #queryseq = queryseq[rtrim:-ltrim]
            #quality_score = quality_score[rtrim:-ltrim]
            # elif rtrim :
            #refseq = refseq[rtrim:]
            #queryseq = queryseq[rtrim:]
            #quality_score = quality_score[rtrim:]

            # remember it is reversed!
            #chromosome_name = self.samfile.getrname(aread.tid)
            for r, q, s in reversed(zip(refseq, queryseq, quality_score)):
                if offset >= counting_start:
                    if offset >= counting_end:
                        break
                    else:
                        if q.isalpha() and s >= cutoff:
                            if r == 'G':
                                methylcontext = self.get_methylation_context(
                                                    roffset, direction=2)
                                if methylcontext == 'CHH':
                                    if q == 'G':
                                        chhm[offset] += 1
                                    elif q == 'A':
                                        chhu[offset] += 1
                                    else:
                                        chho[offset] += 1
                                elif methylcontext == 'CHG':
                                    if q == 'G':
                                        chgm[offset] += 1
                                    elif q == 'A':
                                        chgu[offset] += 1
                                    else:
                                        chgo[offset] += 1
                                else:  # methylcontext == 'CpG' :
                                    if q == 'G':
                                        cpgm[offset] += 1
                                    elif q == 'A':
                                        cpgu[offset] += 1
                                    else:
                                        cpgo[offset] += 1

                if q.isalpha():
                    offset += 1
                if r.isalpha():
                    roffset -= 1
                    assert(roffset >= 0)

    def count_methylation_directional_single(self, aread):
        rtrim = self.rtrim
        ltrim = self.ltrim
        
        #determine the left (5' end) and right (3' end) trims given in nucleotides
        #as in the index.
        if aread.is_reverse :
            ltrim, rtrim = rtrim, ltrim
        
        trim_start, trim_end = self.determine_trimmed_range(aread)
        
        counting_start = max(0, ltrim, trim_start)
        counting_end = min( aread.rlen, aread.rlen - rtrim, trim_end )    
        
        self._count_methylation_simple(aread, counting_start, counting_end)
        
        
    def _count_methylation_simple(self, aread, iqual, start=None, end=None, strand=None):
        '''
        interval method for final stage of the counting.
        aread is the AlignedRead object from pysam.
        start and end can be defined for incorporating the trimming methods that 
        determines the boundaries of good portion or read either from 5' end or
        from 3' end. More delicate nucleotide counting exclusion by quality score
        can be done by class property "quality_score_cutoff".
        
        Note that the strand option is there for dealing with directional paired end
        sequencing data. The paired end sequencing reads can contain the information
        of the strand that is not the same as the mapped strand of the reads.
        '''
        
        if self.verbose > 2 :
            print >>stderr, '*'*20
            print >>stderr, aread.qname, aread.pos, aread.aligned_pairs
            print >>stderr, aread.seq, aread.is_reverse
            print >>stderr, strand
            print >>stderr, '*'*20
            
        if not strand :
            strand = '+'
            if aread.is_reverse :
                strand = '-'
        
        counts = self.counts[strand]
        a, c, g, t, n = counts['A'], counts['C'], counts['G'], counts['T'], counts['N']
        cutoff = self.quality_score_cutoff
        
        for (i,r), s, q in zip(aread.aligned_pairs[start:end], 
                             aread.seq[start:end], 
                             iqual[start:end]) :
            
            #cannot count if the position is unaligned
            #or insertion compared to the reference
            if i == None or r == None :
                continue
            
            #if the position is bad quality 
            #according to the exclusion quality 
            if q <= cutoff :
                continue
            
            if s == 'A' or s == 'a' : a[r] += 1
            elif s == 'T' or s == 't' : t[r] += 1
            elif s == 'G' or s == 'g' : g[r] += 1
            elif s == 'C' or s == 'c' : c[r] += 1
            else : n[r] += 1
        
            if self.verbose > 2 :
                print >>stderr, 'i:', i, 'r:', r, 's:', s, 'q:', q
                print >>stderr, 'ATGC:', a[r], t[r], g[r], c[r]
             
        
    def count_methylation(self, aread, refseq, queryseq, quality_score):
        # update the class variables using the aligned refseq and queryseq
        # only counts with upper characters!! since soft clipped regions supposed
        # to be lower characters.
        rtrim = self.rtrim
        ltrim = self.ltrim
        cutoff = self.quality_score_cutoff

        # new way of incorporating ltrim and rtrim (5' and 3' end trimming)
        counting_start = ltrim
        counting_end = aread.rlen - rtrim

        if self.quality_trim != None:
            counting_end = min(
                counting_end, self.determine_trim_index(aread, quality_score, queryseq))

        # forward strand query
        if not aread.is_reverse:
            pa, pt, pc, pg = self.counts[
                '+']['A'], self.counts['+']['T'], self.counts['+']['C'], self.counts['+']['G']
            offset = aread.qstart
            roffset = aread.pos

            # if rtrim :
            # cut few nt at the 3' terminal
            #refseq = refseq[ltrim:-rtrim]
            #queryseq = queryseq[ltrim:-rtrim]
            #quality_score = quality_score[ltrim:-rtrim]

            for r, q, s in zip(refseq, queryseq, quality_score):
                if offset >= counting_start:
                    if offset >= counting_end:
                        break
                    else:
                        if q.isalpha() and s >= cutoff:
                            if q == 'T':
                                pt[roffset] += 1
                            elif q == 'C':
                                pc[roffset] += 1
                            elif q == 'G':
                                pg[roffset] += 1
                            elif q == 'A':
                                pa[roffset] += 1

                if q.isalpha():
                    offset += 1

                if r.isalpha():
                    roffset += 1  # offset is the reference index.
        else:
            na, nt, nc, ng = self.counts[
                '-']['A'], self.counts['-']['T'], self.counts['-']['C'], self.counts['-']['G']
            # read start position for the backward case!
            # Note that the read is reversed!
            offset = aread.rlen - aread.qend  # backward adjustment
            roffset = aread.aend - 1  # adjustment for reversing

            assert(offset >= 0)
            # if trim :
            # cut few nt at the 3' terminal
            #refseq = refseq[trim:]
            #queryseq = queryseq[trim:]
            #quality_score = quality_score[trim:]
            #roffset = roffset + trim

            for r, q, s in reversed(zip(refseq, queryseq, quality_score)):
                if offset >= counting_start:
                    if offset >= counting_end:
                        break
                    else:
                        if q.isalpha() and s >= cutoff:
                            if q == 'A':
                                na[roffset] += 1
                            elif q == 'G':
                                ng[roffset] += 1
                            elif q == 'T':
                                nt[roffset] += 1
                            elif q == 'C':
                                nc[roffset] += 1

                if q.isalpha():
                    offset += 1

                if r.isalpha():
                    roffset -= 1

    
    def ri2ci(self, aread, read_range):
        '''
        Converts a raw read index (the 0 based index in the given fastq record)
        into a chromosome index (the 0 based index based on the reference chromosome)
        '''
        s, e = read_range
        for sqi, sci in aread.aligned_pairs[s:e] :  # @UnusedVariable
            if sci :
                break
        else :
            raise ChromosomeIndexOutOfBoundError(aread, read_range, aread.aligned_pairs[s:e])
        
        for eqi, eci in reversed(aread.aligned_pairs[s:e]) :  # @UnusedVariable
            if eci :
                break
        else :
            raise ChromosomeIndexOutOfBoundError(aread, read_range, aread.aligned_pairs[s:e])
        
        return sci, eci+1
    
    
    def ci2ri(self, aread, chr_range):
        '''
        Reverse conversion of the ri2ci method.
        '''
        sci, eci = chr_range
        for ri, (qi, ci) in enumerate(aread.aligned_pairs) :  # @UnusedVariable
            if ci >= sci :
                sri = ri
                break
        else :
            #raise RawReadIndexOutOfBoundError(aread)
            sri = ri
        
        for ri, (qi, ci) in reversed(list(enumerate(aread.aligned_pairs))) :  # @UnusedVariable
            if ci < eci :
                eri = ri
                break
        else :
            #raise RawReadIndexOutOfBoundError(aread)
            eri = ri

        if sri > eri :
            return sri, sri
        else :
            return sri, eri + 1

    
    
    def determine_mid_point(self, aread1, p1, aread2, p2):
        
        #take care of singular cases of 0 length indices
        if p1[0] == p1[1] and p2[0]==p2[1] :
            return (0,0), (0,0)
        elif p1[0] == p1[1] :
            return (0,0), p2
        elif p2[0] == p2[1] :
            return p1, (0,0)
        
        cp1 = self.ri2ci(aread1, p1)
        cp2 = self.ri2ci(aread2, p2)
        try :
            if overlap( cp1, cp2 ) :
                #supposedly find if the two ranges are overlapping
                #but if one is bigger than the other, the function should raise an
                #exception
                #the two ranges cp1, cp2 can be exactly same
                cp1, cp2 = split_middle( cp1, cp2 )
                return self.ci2ri(aread1, cp1), self.ci2ri(aread2, cp2)
            else :
                return p1, p2
        except RangeInclusionError :
            if p1[0] <= p2[0] and p1[1] >= p2[1] :
                return p1, (0,0)
            else :
                return (0,0), p2
        
    def count_methylation_directional_pair(self, aread1, aread2, iqual1, iqual2):
        rtrim = self.rtrim
        ltrim = self.ltrim
        
        #determine the left (5' end) and right (3' end) trims given in nucleotides
        #as in the index.
        if aread1.is_reverse :
            ltrim1, rtrim1 = rtrim, ltrim
        else :
            rtrim1, ltrim1 = rtrim, ltrim
            
        if aread2.is_reverse :
            ltrim2, rtrim2 = rtrim, ltrim
        else :
            rtrim2, ltrim2 = rtrim, ltrim
        
        trim_start1, trim_end1 = self.determine_trimmed_range(aread1, iqual1)
        counting_start1 = max(0, ltrim1, trim_start1)
        counting_end1 = min( aread1.rlen, aread1.rlen - rtrim1, trim_end1 )    
        
        trim_start2, trim_end2 = self.determine_trimmed_range(aread2, iqual2)
        counting_start2 = max(0, ltrim2, trim_start2)
        counting_end2 = min( aread2.rlen, aread2.rlen - rtrim2, trim_end2 )    
        
        #additional step in paired end reads
        (counting_start1, counting_end1), (counting_start2, counting_end2) = \
            self.determine_mid_point(aread1, (counting_start1, counting_end1), 
                                     aread2, (counting_start2, counting_end2))
        
        #strand information is necessary for the directional paired end reads
        #has information for the first read.
        if aread1.is_reverse :
            strand = '-'
        else :
            strand = '+'
            
        self._count_methylation_simple(aread1, iqual1, counting_start1, counting_end1, strand=strand)
        self._count_methylation_simple(aread2, iqual2, counting_start2, counting_end2, strand=strand)

                            
    def get_methyltarget(self, i):
        major_cutoff = self.snp_major_allele_cutoff
        min_depth = self.snp_detection_min_depth
        # helper function to determine the genotype
        gp = self.counts['+']['G'][i]
        gn = self.counts['-']['C'][i]
        np = self.counts['+']['A'][i] + \
            self.counts['+']['C'][i] + gp + self.counts['+']['T'][i]
        nn = self.counts['-']['A'][i] + gn + \
            self.counts['-']['G'][i] + self.counts['-']['T'][i]

        if np >= min_depth and gp * 1.0 / np > major_cutoff:
            return '-'
        elif nn >= min_depth and gn * 1.0 / nn > major_cutoff:
            return '+'
        else:
            return ''


    def bsmap_correction(self, totalC, i, strand):
        # need to multiply the factor to the observed C+T
        # i: position in the chromosome
        if strand == '+' or strand == 1:
            a = self.counts['-']['T'][i]
            g = self.counts['-']['C'][i]
        else:
            a = self.counts['+']['A'][i]
            g = self.counts['+']['G'][i]

        try:
            return totalC * g * 1.0 / (a + g)
        except:
            return totalC


    def allnt_correction(self, totalC, i, strand):
        # need to be subtracted from the observed C+T
        # i: position in the chromosome
        if strand == '+' or strand == 1:
            np = self.counts['+']['A'][i] + self.counts['+']['C'][i] + \
                self.counts['+']['G'][i] + self.counts['+']['T'][i]
            an = self.counts['-']['T'][i]
            nn = an + \
                self.counts['-']['C'][i] + \
                self.counts['-']['G'][i] + self.counts['-']['A'][i]
        else:
            np = self.counts['-']['A'][i] + self.counts['-']['C'][i] + \
                self.counts['-']['G'][i] + self.counts['-']['T'][i]
            an = self.counts['+']['A'][i]
            nn = an + \
                self.counts['+']['C'][i] + \
                self.counts['+']['G'][i] + self.counts['+']['T'][i]

        try:
            return totalC - an * 1.0 / nn * np
        except:
            return totalC

    def build_mbiasplot_output(self, sfp=None):
        if not sfp:
            sfp = StringIO()

        #pm, pu, po = self.readcounts['+']['M'], self.readcounts['+']['U'], self.readcounts['+']['O']
        cpgpm, cpgpu, cpgpo = self.readcounts_contexts['CpG'][
            '+']['M'], self.readcounts_contexts['CpG']['+']['U'], self.readcounts_contexts['CpG']['+']['O']
        chgpm, chgpu, chgpo = self.readcounts_contexts['CHG'][
            '+']['M'], self.readcounts_contexts['CHG']['+']['U'], self.readcounts_contexts['CHG']['+']['O']
        chhpm, chhpu, chhpo = self.readcounts_contexts['CHH'][
            '+']['M'], self.readcounts_contexts['CHH']['+']['U'], self.readcounts_contexts['CHH']['+']['O']

        cpgnm, cpgnu, cpgno = self.readcounts_contexts['CpG'][
            '-']['M'], self.readcounts_contexts['CpG']['-']['U'], self.readcounts_contexts['CpG']['-']['O']
        chgnm, chgnu, chgno = self.readcounts_contexts['CHG'][
            '-']['M'], self.readcounts_contexts['CHG']['-']['U'], self.readcounts_contexts['CHG']['-']['O']
        chhnm, chhnu, chhno = self.readcounts_contexts['CHH'][
            '-']['M'], self.readcounts_contexts['CHH']['-']['U'], self.readcounts_contexts['CHH']['-']['O']
        #nm, nu, no = self.readcounts['-']['M'], self.readcounts['-']['U'], self.readcounts['-']['O']

        print >>sfp, "Position\tCpGmC+\tCpGmC-\tCpGC+\tCpGC-\tCpGE+\tCpGE-\tCpGMR+\tCpGMR-\tCpGER+\tCpGER-\tCHGmC+\tCHGmC-\tCHGC+\tCHGC-\tCHGE+\tCHGE-\tCHGMR+\tCHGMR-\tCHGER+\tCHGER-\tCHHmC+\tCHHmC-\tCHHC+\tCHHC-\tCHHE+\tCHHE-\tCHHMR+\tCHHMR-\tCHHER+\tCHHER-"
        for i, (cpgm1, cpgm2, cpgu1, cpgu2, cpgo1, cpgo2, chgm1, chgm2, chgu1, chgu2, chgo1, chgo2, chhm1, chhm2, chhu1, chhu2, chho1, chho2) in enumerate(zip(cpgpm, cpgnm, cpgpu, cpgnu, cpgpo, cpgno,  chgpm, chgnm, chgpu, chgnu, chgpo, chgno,  chhpm, chhnm, chhpu, chhnu, chhpo, chhno)):
            if cpgm1 + cpgm2 + cpgu1 + cpgu2 + cpgo1 + cpgo2 + chgm1 + chgm2 + chgu1 + chgu2 + chgo1 + chgo2 + chhm1 + chhm2 + chhu1 + chhu2 + chho1 + chho2:
                try:
                    cpgt1, cpgt2, chgt1, chgt2, chht1, chht2 = 1.0 * cpgm1 + cpgu1, 1.0 * cpgm2 + \
                        cpgu2, 1.0 * chgm1 + chgu1, 1.0 * chgm2 + \
                        chgu2, 1.0 * chhm1 + chhu1, 1.0 * chhm2 + chhu2

                    cpgmr1, cpgmr2, chgmr1, chgmr2, chhmr1, chhmr2 = cpgm1 / cpgt1, cpgm2 / \
                        cpgt2, chgm1 / chgt1, chgm2 / \
                        chgt2, chhm1 / chht1, chhm2 / chht2
                    #cpgt1, cpgt2, chgt1, chgt2, chht1, chht2 = cpgt1+cpgo1, cpgt2+cpgo2, chgt1+chgo1, chgt2+chgo2, chht1+chho1, chht2+chho2
                    cpgor1, cpgor2, chgor1, chgor2, chhor1, chhor2 = cpgo1 / cpgt1, cpgo2 / \
                        cpgt2, chgo1 / chgt1, chgo2 / \
                        chgt2, chho1 / chht1, chho2 / chht2
                except ZeroDivisionError:
                    print >>sfp, '\t'.join([str(v) for v in (i + 1, cpgm1, cpgm2, cpgu1, cpgu2, cpgo1, cpgo2, '-', '-', '-', '-', chgm1,
                                                             chgm2, chgu1, chgu2, chgo1, chgo2, '-', '-', '-', '-',  chhm1, chhm2, chhu1, chhu2, chho1, chho2, '-', '-', '-', '-',)])
                    continue

                print >>sfp, '\t'.join([str(v) for v in (i + 1, cpgm1, cpgm2, cpgu1, cpgu2, cpgo1, cpgo2, cpgmr1, cpgmr2, cpgor1, cpgor2, chgm1, chgm2,
                                                         chgu1, chgu2, chgo1, chgo2, chgmr1, chgmr2, chgor1, chgor2, chhm1, chhm2, chhu1, chhu2, chho1, chho2, chhmr1, chhmr2, chhor1, chhor2)])

    def build_methylation_bedgraph(self, sfp=None, snp=0, comprehensive_output=True):
        # sfp: file pointer to print out the output
        # snp: snp correction type flag,
            # 0-No correction,
            # 1-All nucleotide correction using A-/(N-)*(N+),
            # 2-BSMAP correction, using G-/(A-+G-)

        if not sfp:
            sfp = StringIO()

        pa, pt, pc, pg = self.counts[  # @UnusedVariable
            '+']['A'], self.counts['+']['T'], self.counts['+']['C'], self.counts['+']['G']
        na, nt, nc, ng = self.counts[  # @UnusedVariable
            '-']['A'], self.counts['-']['T'], self.counts['-']['C'], self.counts['-']['G']

        if snp:
            for i, r in self.iterposnt():
                # r: nucleotide in the reference genome
                # i: 0-based position in the chromosome
                methyltarget = self.get_methyltarget(i)

                effective_totalC = 0.0
                effective_ncount = 0.0

                if methyltarget == '+':
                    mcount = pc[i]
                    ncount = pt[i]
                    totalC = 0.0 + mcount + ncount
                    genotype = 'C'
                elif methyltarget == '-':
                    mcount = ng[i]
                    ncount = na[i]
                    totalC = 0.0 + mcount + ncount
                    genotype = 'G'
                elif r == 'c' or r == 'C':
                    mcount = pc[i]
                    ncount = pt[i]
                    totalC = 0.0 + mcount + ncount
                    genotype = 'N'
                elif r == 'g' or r == 'G':
                    mcount = ng[i]
                    ncount = na[i]
                    totalC = 0.0 + mcount + ncount
                    genotype = 'N'
                else:
                    continue  # skip unless the position is methyl target

                if totalC > 0.0:
                    if snp == 1:
                        effective_totalC = self.bsmap_correction(
                            totalC, i, methyltarget)
                    elif snp == 2:
                        effective_totalC = self.allnt_correction(
                            totalC, i, methyltarget)
                    else:
                        effective_totalC = totalC

                effective_ncount = max(effective_totalC - mcount, 0.0)

                if not comprehensive_output:
                    print >>sfp, '\t'.join(
                        [self.name, str(i), str(i + 1), str(mcount), str(effective_ncount)])
                else:
                    print >>sfp, '\t'.join([self.name, str(i), str(
                        i + 1), str(mcount), str(effective_ncount), str(ncount), genotype, r])

        else:
            for i, r in self.itermethylposseq():
                # r: nucleotide in the reference genome
                # i: 0-based position in the chromosome
                if r == 'c' or r == 'C':
                    print >>sfp, '\t'.join(
                        [self.name, str(i), str(i + 1), str(pc[i]), str(pt[i])])
                else:
                    print >>sfp, '\t'.join(
                        [self.name, str(i), str(i + 1), str(ng[i]), str(na[i])])
        return sfp


    def build_methylation_summary(self, sfp=None):
        # sfp: file pointer to print out the output
        # snp: snp correction type flag,
            # 0-No correction,
            # 1-All nucleotide correction using A-/(N-)*(N+),
            # 2-BSMAP correction, using G-/(A-+G-)

        if not sfp:
            sfp = StringIO()

        pa, pt, pc, pg = self.counts[  # @UnusedVariable
            '+']['A'], self.counts['+']['T'], self.counts['+']['C'], self.counts['+']['G']
        na, nt, nc, ng = self.counts[  # @UnusedVariable
            '-']['A'], self.counts['-']['T'], self.counts['-']['C'], self.counts['-']['G']

        #total counts
        tc = {'+':{'CpG':0, 'CHG':0, 'CHH':0}, '-':{'CpG':0, 'CHG':0, 'CHH':0}}
        #methylation counts
        mc = {'+':{'CpG':0, 'CHG':0, 'CHH':0}, '-':{'CpG':0, 'CHG':0, 'CHH':0}}

        for i, direction, context in self.itermethylitems():
            # r: nucleotide in the reference genome
            # i: 0-based position in the chromosome
            if direction == '+':
                #print >>sfp, '\t'.join(
                #    [self.name, str(i), str(i + 1), str(pc[i]), str(pt[i])])
                tc[direction][context] += pc[i]+pt[i]
                mc[direction][context] += pc[i]
            else:
                #print >>sfp, '\t'.join(
                #    [self.name, str(i), str(i + 1), str(ng[i]), str(na[i])])
                tc[direction][context] += na[i]+ng[i]
                mc[direction][context] += ng[i]
        
        
        for direction in ('+', '-') :
            for context in ('CpG', 'CHG', 'CHH') :
                print >>sfp, '\t'.join([self.name, str(mc), str(tc)])
                
        return sfp
