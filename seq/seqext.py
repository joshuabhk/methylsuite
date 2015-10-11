import cStringIO
from random import getrandbits
from seq import reverse_complement

class SequenceExtracter :
        def __init__( self, fa ) :
                #setting up FASTA file for extracting chromosomes
                self.chromosomes = {}
                if fa :
                        fp = open( fa )
                        strfp = cStringIO.StringIO()
                        for l in fp :
                                if l.startswith( '>' ) :
                                        #if strfp contains chromosomes
                                        if strfp.tell() :
                                                self.chromosomes[seqid] = strfp.getvalue()  # @UndefinedVariable
                                                strfp = cStringIO.StringIO()
                                        seqid = l[1:].strip()
                                else :
                                        strfp.write( l.strip() )
                        else :
                                if strfp.tell() :
                                        self.chromosomes[seqid] = strfp.getvalue()
                                        strfp = None


        def get_sequence( self, chrom, start, length, direction=None ) :
                #direction can be 0, 1, 2.
                #0 means random selection
                #1 means positive strand
                #2 means negative strand
                if not direction :
                        direction = getrandbits(1)+1

                if direction == 1 or direction == '+' :
                        return self.chromosomes[chrom][start:start+length], '+'
                else :
                        return reverse_complement( self.chromosomes[chrom][start:start+length] ), '-'


