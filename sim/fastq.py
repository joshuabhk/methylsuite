import sys
import cStringIO
from seq import FASTQ
import random

class FASTQRandomizer( FASTQ ) :
        def build_random_quality_string_based_on_input_fastq( self, length, base ) :
                #build sanger range score

                scores = []
                for i in xrange( length ) :
                        if i < len(self.means) :
                                score = self.force_sanger_range(gauss( self.means[i], self.stdevs[i] )) #need to check
                        else :
                                score = self.force_sanger_range(gauss( self.means[-1], self.vars[-1] ))
                        scores.append( chr(score + base) )

                return ''.join( scores )

        def build_randomly_methylated_sequence( self, seq, methylation_rate ) :
                strfp = cStringIO.StringIO()
                for a in seq :
                        if a == 'c' :
                                if random.random() > methylation_rate :
                                        strfp.write( a )
                                else :
                                        strfp.write( 't' )
                        elif a == 'C' :
                                if random.random() > methylation_rate :
                                        strfp.write( a )
                                else :
                                        strfp.write( 'T' )
                        else :
                                strfp.write( a )
                return strfp.getvalue()

        def build_randomly_methylated_records( self, methylation_rate, strfp=None ) :
                if not strfp :
                        stfp = cStringIO.StringIO()

                rec = self.get_record()
                while rec :
                        h1, seq, h2, qual = rec
                        seq = self.build_randomly_methylated_sequence( seq, methylation_rate )

                        print >>strfp, h1
                        print >>strfp, seq
                        print >>strfp, h2
                        print >>strfp, qual

                        rec = self.get_record()

                return strfp


        def build_position_based_methylated_sequence( self, seq, methylrates, offset, strand ) :
                strfp = cStringIO.StringIO()
                for i, a in enumerate( seq ) :
                        if a == 'c' or a == 'C' :
                                if strand == '+' or strand == 1 :
                                        linpos = offset + i
                                else :
                                        linpos = offset + len(seq) - 1 - i

                                try :
                                        methylation_rate = methylrates.linpos2methylrate( linpos )
                                except KeyError :
                                        #SNP Cytosine #into three groups
                                        if linpos % 3 == 0 : #methylation rate 0
                                                methylation_rate = 0.05
                                        elif linpos % 3 == 1 : #  0.5
                                                methylation_rate = 0.05

                                        else : #even number methylation rate 0
                                                methylation_rate = 0.05

                                if a.islower() :
                                        if random.random() < methylation_rate :
                                                strfp.write( a )
                                        else :
                                                strfp.write( 't' )
                                else :
                                        if random.random() < methylation_rate :
                                                strfp.write( a )
                                        else :
                                                strfp.write( 'T' )
                        else :
                                strfp.write( a )

                return strfp.getvalue()

        def build_position_based_methylated_records( self, methylrates, strfp=None ) :
                if not strfp :
                        stfp = cStringIO.StringIO()

                rec = self.get_record()
                while rec :
                        h1, seq, h2, qual = rec
                        readid, chr, pos, length, strand =  parse_header(h1)
                        offset = methylrates.chromosomepos2linpos( chr, int(pos) )

                        seq = self.build_position_based_methylated_sequence( seq, methylrates, offset, strand )

                        print >>strfp, h1
                        print >>strfp, seq
                        print >>strfp, h2
                        print >>strfp, qual

                        rec = self.get_record()

                return strfp


        def introduce_random_methylation_per_nucleotide( self, nucleotide='C', rate=0.5 ) :
                if random.random() > rate :
                        if nucleotide.isupper() :
                                return 'T'
                        else :
                                return 't'
                else :
                        return nucleotide

        def get_cytosine_context( self, strand, chr, pos, seq, genome=None ) :
                seq = seq.upper()
                if len(seq) == 3 :
                        if seq.startswith( 'CG' ) :
                                return 'CpG'
                        elif chg_regex.match( seq ) :
                                return 'CHG'
                        elif chh_regex.match( seq ) :
                                return 'CHH'
                        else :
                                raise Exception( "context could not be determined", seq )

                elif len(seq) == 2 :
                        if seq.startswith( 'CG' ) :
                                return 'CpG'
                        else :
                                if genome :
                                        return genome.get_cytosine_context( strand, chr, pos )
                                else :
                                        raise Exception( "Cannot determine the context" )
                else :
                        if genome :
                                return genome.get_cytosine_context( strand, chr, pos )
                        else :
                                raise Exception( "Cannot determine the context" )


        def build_context_dependent_methylated_sequence( self, strand, chr, pos, sequence, genome=None, rates = {'CpG':0.8, 'CHG':0.1, 'CHH':0.05} ) :
                newseq = []
                for i, a in enumerate( sequence ) :
                        if a == 'C' or a == 'c' :
                                if direction == '+' or direction == 1 :
                                        ccontext = self.get_cytosine_context( strand, chr, pos+i, sequence[i:i+3], genome=genome )
                                else :
                                        ccontext = self.get_cytosine_context( strand, chr, pos+len(sequence)-i-1, sequence[i:i+3], genome=genome )
                                newseq.append( self.introduce_random_methylation_per_nucleotide( nucleotide=a, rate=rates[ccontext] ) )
                        else :
                                newseq.append( a )

                return ''.join( newseq )

        def build_methylated_records( self, genome=None, rates = {'CpG':0.8, 'CHG':0.1, 'CHH':0.05}, strfp=None ) :
                if not strfp :
                        strfp = cStringIO.StringIO()

                rec = self.get_record()
                while rec :
                        h1, seq, h2, qual = rec
                        serial, chr, pos, length, strand = parse_header(h1)

                        seq = self.build_context_dependent_methylated_sequence( strand, chr, int(pos), seq, genome=genome, rates=rates )
                        print >>strfp, h1
                        print >>strfp, seq
                        print >>strfp, h2
                        print >>strfp, qual

                        rec = self.get_record()

                return strfp

        def introduce_random_sequencing_error( self, seq, error_rate=None, quality_scores=None ) :
                newseq = []
                if error_rate != None :
                        for a in seq :
                                if a == 'N' or a == 'n' :
                                        newseq.append(a)

                                elif random.random() < error_rate :
                                        newseq.append( random_nucleotide(a) )
                                else :
                                        newseq.append(a)
                        return ''.join( newseq )

                elif quality_scores :
                        for i, a in enumerate(seq) :
                                if a == 'N' or a == 'n' :
                                        newseq.append(a)
                                elif random.random() < self.quality_score2error_rate( quality_scores[i] ) :
                                        newseq.append( random_nucleotide(a) )
                                else :
                                        newseq.append(a)
                        return ''.join( newseq )

                else :
                        #if nothing given
                        #return as is
                        return seq

        def introduce_random_sequencing_error_into_records( self, error_rate=None, strfp=None ) :
                #if error_rate ([0.0, 1.0] is given, sequencing errors will be introduced with equal rate
                #regardless of the quality score.

                if not strfp :
                        strfp = cStringIO.StringIO()

                rec = self.get_record()
                while rec :
                        h1, seq, h2, qual = rec
                        if error_rate != None :
                                seq = self.introduce_random_sequencing_error( seq, error_rate=error_rate )
                        else :
                                seq = self.introduce_random_sequencing_error( seq, quality_scores=self.quality_string2score(qual) )

                        print >>strfp, h1
                        print >>strfp, seq
                        print >>strfp, h2
                        print >>strfp, qual

                        rec = self.get_record()

                return strfp

        def introduce_random_sequencing_error2( self, seq, error_rate=None, quality_scores=None ) :
                newseq = []
                if error_rate != None :
                        for a in seq :
                                if a == 'N' or a == 'n' :
                                        newseq.append(a)

                                elif random.random() < error_rate :
                                        newseq.append( random_nucleotide(a) )
                                else :
                                        newseq.append(a)
                        return ''.join( newseq )

                elif quality_scores :
                        for i, a in enumerate(seq) :
                                if a == 'N' or a == 'n' :
                                        newseq.append(a)
                                elif choice_with_N( quality_scores[i] ) == 'N' :
                                        if a.isupper() :
                                                newseq.append( 'N' )
                                        else :
                                                newseq.append( 'n' )
                                elif random.random() < self.quality_score2error_rate( quality_scores[i] ) :
                                        newseq.append( random_nucleotide2(a) )
                                else :
                                        newseq.append(a)
                        return ''.join( newseq )

                else :
                        #if nothing given
                        #return as is
                        return seq

        def introduce_random_sequencing_error_into_records2( self, error_rate=None, strfp=None ) :
                #this is second version of the sequencing error introducer.
                #This is improved in modeling errors N and more realistic error models.
                #if error_rate ([0.0, 1.0] is given, sequencing errors will be introduced with equal rate
                #regardless of the quality score.

                if not strfp :
                        strfp = cStringIO.StringIO()

                for rec in self.iterrecords() :
                        h1, seq, h2, qual = rec
                        if error_rate != None :
                                seq = self.introduce_random_sequencing_error( seq, error_rate=error_rate )
                        else :
                                seq = self.introduce_random_sequencing_error2( seq, quality_scores=self.quality_string2score(qual) )

                        print >>strfp, h1
                        print >>strfp, seq
                        print >>strfp, h2
                        print >>strfp, qual

                return strfp

        def introduce_snp( self, snp, chr, pos, seq, strand ) :
                new_seq = []
                for offset, a in enumerate( seq ) :
                        try :
                                new_seq.append( snp.get_random_allele( chr, pos+offset, a, strand ) )
                        except :
                                print >>sys.stderr, chr, pos, offset, a, strand, "has a problem."
                                new_seq.append( a )

                return ''.join( new_seq )

        def introduce_snp_into_records( self, snp, strfp=None ) :
                #snp should be a DBSNPRandomizer object
                #note that this function depends on the chromosome location information
                #in the header line
                #ex) ^<id> <chr> <start> <length> <strand>$
                #ex) ^12 chr1 123421 82 +$

                if not strfp :
                        strfp = cStringIO.StringIO()

                rec = self.get_record()
                while rec :
                        h1, seq, h2, qual = rec

                        #parsing header
                        readid, chr, pos, length, strand = parse_header(h1)

                        seq = self.introduce_snp( snp, chr, int(pos), seq, strand )
                        print >>strfp, h1
                        print >>strfp, seq
                        print >>strfp, h2
                        print >>strfp, qual

                        rec = self.get_record()

                return strfp
