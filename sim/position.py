from random import randrange
from seq.seqext import SequenceExtracter
import cStringIO
from seq.fastq import FASTQ

class DNAPositionSelector :
        def __init__( self, chrominfo, subset=None ) :
                fp = open( chrominfo )
                self.chromosomes = []
                self.lengths = []
                self.ends = [] #end of the chromosomes in the linear position
                for l in fp :
                        chr, length = l.split()  # @ReservedAssignment

                        if subset and not chr in subset :
                                continue

                        self.chromosomes.append( chr )
                        self.lengths.append( int(length) )

                for i, length in enumerate( self.lengths ) :
                        if i -1 > 0 :
                                self.ends.append( self.ends[i-1] + self.lengths[i] )
                        else :
                                #beginning
                                self.ends.append( self.lengths[i] )

        def chromosome2length( self, chr ) :  # @ReservedAssignment
                return self.lengths[ self.chromosomes.index(chr) ]

        def linpos2chromosomepos( self, linpos ) :
                #convert a linear position into a chromosome and position in the chromosome
                for i, end in enumerate( self.ends ):
                        if linpos < end :
                                if i == 0 :
                                        offset = 0
                                else :
                                        offset = self.ends[i-1]

                                return self.chromosomes[i], linpos - offset

                        raise Exception( "No position found" )

        def chromosomepos2linpos( self, chromosome, pos ) :
                i = self.chromosomes.index( chromosome )
                if i == 0 :
                        offset = 0
                else :
                        offset = self.ends[i-1]

                return offset + int(pos)

        def distance_to_end( self, chr, pos ) :  # @ReservedAssignment
                return self.chromosome2length(chr) - pos -1


class PositionSpecificMethylationRates( DNAPositionSelector ) :
        def __init__( self, chrominfo="/home/kimb/projects/comparing_bisulfite_reads/hg19.chromInfo", bedgraph="/home/kimb/isilon_temp/methylation_calls/roadmap/h1nueral/bismark/chr21/chr21.bedGraph" ) :
                #read bedgraph file containing methylation rate information.
                self.methylation_rates = {}
                self.percentage_input = 1
                if chrominfo :
                        DNAPositionSelector.__init__( self, chrominfo=chrominfo )
                if bedgraph :
                        self._set_methylation_rates( bedgraph )

        def _set_methylation_rates( self, bedgraph ) :
                fp = open( bedgraph )
                convert_percentage = self.percentage_input
                for l in fp :
                        line = l.split('\t')
                        chr = line[0]  # @ReservedAssignment
                        start = int(line[1])
                        #end = int(line[2])

                        #Wrong!
                        #bismark adjustment
                        #if start == end :
                        #       start -= 1

                        val = float(line[3])
                        if convert_percentage : #important.
                                val = val/100

                        self.methylation_rates[ self.chromosomepos2linpos( chr, start ) ] = val

        def linpos2methylrate( self, linpos ) :
                return self.methylation_rates[ linpos ]

        def chrpos2methylrate( self, chr, pos ) :  # @ReservedAssignment
                return self.methylation_rates[ self.chrmosomepos2linpos(chr,pos) ]

class RandomSequenceGenerator(DNAPositionSelector) :
        def __init__( self, chrominfo, fasta, subset=[] ) :

                DNAPositionSelector.__init__( self, chrominfo, subset=subset )
                self.sequence_extracter = SequenceExtracter( fasta )
                self.get_sequence = self.sequence_extracter.get_sequence

        def get_random_position(self) :
                linpos = randrange( self.ends[-1] )
                return self.linpos2chromosomepos( linpos )

        def get_random_sequence( self, length=82, direction=None ) :
                chr, pos = self.get_random_position()  # @ReservedAssignment
                #prevent the sequence index out of bound.
                while self.distance_to_end( chr, pos ) < length :
                        chr, pos = self.get_random_position()  # @ReservedAssignment

                return self.get_sequence( chr, pos, length, direction )
        def get_random_fasta( self, seqid=None, length=82, direction=None, sep='_' ) :
                chr, pos = self.get_random_position()  # @ReservedAssignment

                #prevent the sequence index out of bound.
                while self.distance_to_end( chr, pos ) < length :
                        chr, pos = self.get_random_position()  # @ReservedAssignment

                seq, direction = self.get_sequence( chr, pos, length, direction )
                if not seqid :
                        seqid = sep.join([ chr, str(pos), str(length), direction] )

                if not seqid.startswith('>') :
                        seqid = '>' + seqid

                return seqid + '\n' + seq + '\n'

        def get_random_fastqs( self, length=82, direction=None, counts=1000000, quality_score_template="", strfp=None, sep='_' ) :
                if not strfp :
                        strfp = cStringIO.StringIO()

                fastq = None
                if quality_score_template : #a usual fastq file whose quality score line will be copied
                        fastq = FASTQ( quality_score_template )
                else :
                        quality_string = FASTQ().build_high_quality_string()

                for i in xrange(counts) :
                        chr, pos = self.get_random_position()  # @ReservedAssignment
                        while self.distance_to_end( chr, pos ) < length :
                                chr, pos = self.get_random_position()  # @ReservedAssignment

                        seq, setdirection = self.get_sequence( chr, pos, length, direction )
                        seqid = sep.join( [str(i), chr, str(pos), str(length), setdirection] )
                        header1 = '@'+sep + seqid
                        header2 = '+'+sep + seqid
                        if fastq :
                                quality_string = fastq.get_next_quality_score_line()

                        print >>strfp, header1
                        print >>strfp, seq
                        print >>strfp, header2
                        print >>strfp, quality_string

                return strfp #note that the return value is File IO object!

        def get_random_fastqs_with_copy_N_pattern( self, length=82, direction=None, counts=1000000, quality_score_template="", strfp=None, sep='_' ) :
                if not strfp :
                        strfp = cStringIO.StringIO()

                fastq = None
                if quality_score_template : #a usual fastq file whose quality score line will be copied
                        fastq = FASTQ( quality_score_template )
                else :
                        quality_string = FASTQ().build_high_quality_string()

                for i in xrange(counts) :
                        chr, pos = self.get_random_position()  # @ReservedAssignment
                        while self.distance_to_end( chr, pos ) < length :
                                chr, pos = self.get_random_position()  # @ReservedAssignment

                        seq, setdirection = self.get_sequence( chr, pos, length, direction )
                        seqid = sep.join( [str(i), chr, str(pos), str(length), setdirection] )
                        header1 = '@'+sep+ + seqid
                        header2 = '+'+sep + seqid
                        if fastq :
                                quality_string = fastq.get_next_quality_score_line()

                        print >>strfp, header1
                        print >>strfp, seq
                        print >>strfp, header2
                        print >>strfp, quality_string

                return strfp #note that the return value is File IO object!


