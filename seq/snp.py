from sim.position import DNAPositionSelector
import gzip

class DBSNP :
        #simple parser for common SNP txt file downloaded from UCSC genome browser.
        def __init__( self, fn, chrominfo_fn ) :

                if fn.endswith( '.gz' ) :
                        fp = gzip.GzipFile(fn)
                else :
                        fp = open(fn)

                self.pc = DNAPositionSelector(chrominfo_fn) # linear coordinate selector
                self.chromosomepos2linpos = self.pc.chromosomepos2linpos
                self.linpos2chromosomepos = self.pc.linpos2chromosomepos

                self.snp = {}
                for l in fp :
                        line = l.strip().split('\t')
                        chr = line[1]  # @ReservedAssignment
                        start = int(line[2])
                        end = int(line[3])
                        type= line[11]  # @ReservedAssignment
                        weight = int( line[17] ) #or weight in the UCSC schema
                        exception = line[18]  # @UnusedVariable

                        if type != 'single' :
                                continue
                        if end-start != 1 :
                                continue
                        if weight > 1 :
                                continue

                        linpos = self.chromosomepos2linpos( chr, start )
                        strand = line[6]
                        alleles = line[22]
                        frequencies = line[24]

                        self.snp[ linpos ] = '\t'.join( [strand, alleles, frequencies] )

        def parse_snp_info( self, str ) :  # @ReservedAssignment
                #internally saved snp information string parser
                strand, allele, frequencies = str.split('\t')
                return strand, allele.split(',')[:-1], [float(f) for f in frequencies.split(',')[:-1] ]

        def get_snp_info( self, chr, pos ) :  # @ReservedAssignment
                linpos = self.chromosomepos2linpos( chr, pos )

                snpinfo = self.snp.get( linpos )
                if snpinfo :
                        snpinfo = self.parse_snp_info( snpinfo )
                        return snpinfo


