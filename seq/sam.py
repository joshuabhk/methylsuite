import bz2
import sys

class SAM :
        '''
        This is a light weighted SAM file parser.
        '''
        def __init__( self, fn='', fp=None ) :
                if fn :
                        if fn.endswith('.bz2' ) :
                                self.fp = bz2.BZ2File(fn)
                        else :
                                self.fp = open(fn)
                else :
                        self.fp = fp

                self.file = fn
                self.mappings = {}

                self.parse()

        def parse_header( self, header ) :
                '''
                Parses FASTQ Header when the fields are concaternated by underscores.
                '''
                fields = header.split('_')
                if len(fields) == 6 :
                        junk, id, chr, pos, length, direction = fields  # @UnusedVariable @ReservedAssignment
                        return id, chr, pos, length, direction
                elif len(fields) > 6 :
                        return fields[1], '_'.join( fields[2:-3] ), fields[-3], fields[-2], fields[-1]
                else :
                        print >>sys.stderr, header
                        assert 0


        def parse( self, unique=True ) :
                #Need to convert to 0 based mapping position for SAM.
                for l in self.fp :
                        if l.startswith( '@' ) :
                                continue
                        line = l.split()
                        self.mappings[ line[0] ] = line[2] + '/' + str(int(line[3])-1)

        def is_mapped( self, readid ) :
                return readid in self.mappings

        def get_mapping( self, readid ) :
                #0 based mapping results should be returned.
                return self.mappings[ readid ].split('/')

        def get_refmapping( self, readid ) :
                return self.parse_header( readid )[1:3]


class SAMCheck(SAM) :
        def check_records( self ) :
                correct, incorrect = 0, 0
                for readid in self.mappings :
                        chr, pos = self.get_refmapping(readid)  # @ReservedAssignment
                        chrmap, posmap = self.get_mapping(readid)

                        if chrmap == chr and  posmap == pos :
                                correct += 1
                        else :
                                incorrect += 1

                        return correct, incorrect

        def check_records_and_compare_mappings( self, readids, sam, remove_non_regular_chromosome=1 ) :
                #perform checking to the reference and comparing between mapping results.
                is_mapped1 = self.is_mapped
                is_mapped2 = sam.is_mapped

                get_mapping1 = self.get_mapping
                get_mapping2 = sam.get_mapping
                get_refmapping = self.get_refmapping

                ac12 = 0 # both methods agreed on mappings & correct
                ai12 = 0 # both methods agreed on mappings & incorrect

                dc1 = 0 # disagree methods & correct only in method 1
                #di1 = 0 # disagree methods & incorrectly only in method 1 #logically same as dc2
                dc2 = 0
                #di2 = 0 # same as dc1
                di12 = 0 #disagree methods & neither of them are right!

                c1 = 0 #correct & method 1 only mapped (meanining itself)
                i1 = 0 #incorrect & method 1 only mapped
                c2 = 0
                i2 = 0
                u12 = 0 #unmapped by both

                t = 0 #total counts
                for readid in readids :
                        #print readid

                        refchr, refpos = get_refmapping( readid )
                        if remove_non_regular_chromosome and '_' in refchr :
                                continue

                        t += 1
                        if is_mapped1( readid ) and is_mapped2(readid) :
                                chr1, pos1 = get_mapping1(readid)
                                chr2, pos2 = get_mapping2(readid)

                                #agreed
                                if chr1 == chr2 and pos1 == pos2 :
                                        if refchr == chr1 and refpos == pos1 :
                                                ac12 += 1
                                        else :
                                                ai12 += 1
                                #disagreed
                                else :
                                        if chr1 == refchr and pos1 == refpos :
                                                dc1 += 1
                                                #di2 += 1
                                        elif chr2 == refchr and pos2 == refpos :
                                                dc2 += 1
                                                #di1 += 1
                                        else :
                                                di12 += 1
                        elif is_mapped1( readid) :
                                chr1, pos1 = get_mapping1(readid)
                                if refchr == chr1 and refpos == pos1 :
                                        c1 += 1
                                else :
                                        i1 += 1

                        elif is_mapped2( readid ) :
                                chr2, pos2 = get_mapping2(readid)
                                if refchr == chr2 and refpos == pos2 :
                                        c2 += 1
                                else :
                                        i2 += 1

                        else :
                                u12 +=1

                return t, ac12, dc1, dc2, ai12, di12, c1, i1, c2, i2, u12

        def check_records_old( self, fastq=None ) :
                if fastq :
                        raise NotImplemented()
                        #get the records from fastq

                self.fp.seek(0)

                correct = 0
                incorrect = 0
                for l in self.fp :
                        if l.startswith( '@' ) :
                                continue
                        line = l.split()
                        chrmap = line[2]
                        posmap = int(line[3])

                        #unmapped read
                        if chrmap == '*' and posmap==0 :
                                continue

                        if fastq :
                                symbol, id, chr, pos, length, direction = fastq.get_record()[0].split( '_' )  # @ReservedAssignment @UnusedVariable
                        else :
                                #print line[0].split('_')
                                '''
                                #wrong! parsing
                                header = line[0].split('_')
                                if len(header) == 6 :
                                        junk, id, chr, pos, length, direction = header
                                elif len(header) > 6 :
                                        id, chr = header[1:3]
                                        pos = header[-3]
                                        length = header[-2]
                                        direction = header[-1]
                                else :
                                        print header
                                        assert 0
                                '''
                                #correct
                                id, chr, pos, length, direction = self.parse_header( line[0] )  # @ReservedAssignment @UnusedVariable

                        if chrmap == chr and  posmap -1 == int(pos) :
                                correct += 1
                        else :
                                incorrect += 1

                return correct, incorrect

class PASHCheck( SAMCheck ) :
        def parse( self ) :
                #convert 1 based position into 0 based.
                #omit multimappers
                multimappers = set()
                id2chrpos = {}
                for l in self.fp :
                        line = l.split('\t')
                        chrpos = line[0] + '/' + str(int(line[1])-1) #conversion
                        id = line[3]  # @ReservedAssignment
                        if id in multimappers :
                                continue
                        elif id in id2chrpos :
                                del id2chrpos[id]
                                multimappers.add( id )
                        else :
                                id2chrpos[id] = chrpos

                self.mappings = id2chrpos

        def check_records_old( self, fastq=None ) :
                if fastq :
                        raise NotImplemented()

                correct = 0
                incorrect = 0

                #checking only uniquely mapped reads
                multimappers = set()
                id2chrpos = {}
                for l in self.fp :
                        line = l.split('\t')
                        chrpos = line[0] + '/' + line[1]
                        id = line[3]  # @ReservedAssignment
                        if id in multimappers :
                                continue
                        elif id in id2chrpos :
                                del id2chrpos[id]
                                multimappers.add( id )
                        else :
                                id2chrpos[id] = chrpos

                for k, v in id2chrpos.iteritems() :
                        '''
                        if len(header) == 6 :
                                junk, id, chr, pos, length, direction = header
                        elif len(header) > 6 :
                                id, chr = header[1:3]
                                pos = header[-3]
                                length = header[-2]
                                direction = header[-1]
                        else :
                                print header
                                assert 0
                        '''
                        #correct
                        id, chr, pos, length, direction = self.parse_header( k )  # @ReservedAssignment @UnusedVariable

                        chrmap, posmap = v.split('/')

                        if chr == chrmap and int(posmap)-1 == int(pos) :
                                correct += 1
                        else :
                                incorrect += 1

                return correct, incorrect


class GsnapSAMCheck( SAMCheck ) :
        def parse( self, unique=True ) :
                #Need to convert to 0 based mapping position for SAM.
                #Check only unique mappers
                for l in self.fp :
                        if l.startswith( '@' ) :
                                continue

                        line = l.split()

                        chrmap = line[2]
                        posmap = int(line[3])

                        #unmapped read
                        if chrmap == '*' and posmap==0 :
                                continue
                        #skipping if the read alinged more than once hits :
                        if line[13].startswith( 'HI:i:' ) :
                                nhits = int(line[13].split(':')[-1])
                                if nhits > 1 :
                                        continue
                        else :
                                print >>sys.stderr, "Error:", "HI tag is missing."
                                print >>sys.stderr, line
                                sys.exit()

                        self.mappings[ line[0] ] = chrmap + '/' + str(posmap-1) #conversion

        def check_records_old( self, fastq=None ) :
                if fastq :
                        raise NotImplemented()
                        #get the records from fastq

                correct = 0
                incorrect = 0
                for l in self.fp :
                        if l.startswith( '@' ) :
                                continue

                        line = l.split()
                        chrmap = line[2]
                        posmap = int(line[3])

                        #unmapped read
                        if chrmap == '*' and posmap==0 :
                                continue

                        #skipping if the read alinged more than once hits :
                        if line[13].startswith( 'HI:i:' ) :
                                nhits = int(line[13].split(':')[-1])
                                if nhits > 1 :
                                        continue
                        else :
                                print >>sys.stderr, "Error:", "HI tag is missing."
                                print >>sys.stderr, line
                                sys.exit()
                        if fastq :
                                symbol, id, chr, pos, length, direction = fastq.get_record()[0].split( '_' )  # @ReservedAssignment @UnusedVariable
                        else :
                                '''
                                #wrong! but not many reads are affected. :)
                                header = line[0].split('_')
                                if len(header) == 6 :
                                        junk, id, chr, pos, length, direction = header
                                elif len(header) > 6 :
                                        id, chr = header[1:3]
                                        pos = header[-3]
                                        length = header[-2]
                                        direction = header[-1]
                                else :
                                        print header
                                        assert 0
                                '''
                                #correct
                                #temporarilly blocked due to the refactoring.
                                #Presumably this method should not be used.
                                #id, chr, pos, length, direction = self.parse_header( k )

                        if chrmap == chr and  posmap -1 == int(pos) :
                                correct += 1
                        else :
                                incorrect += 1

                return correct, incorrect

