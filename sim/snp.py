from seq.snp import DBSNP
from seq import complement
from random import random

class DBSNPRandomizer( DBSNP ) :
        def get_random_allele( self, chr, pos, refallele, refstrand ) :  # @ReservedAssignment
                snpinfo = self.get_snp_info( chr, pos )

                if snpinfo :
                        strand, alleles, frequencies = snpinfo
                        #debug
                        #print strand, alleles, frequencies
                        r = random()
                        f = [0.0]*len(frequencies)
                        selected_index = len(frequencies)-1

                        if selected_index < 0 :
                                raise Exception( selected_index )

                        for i, freq in enumerate(frequencies) :
                                if i == 0 :
                                        f[i] = freq
                                else :
                                        f[i] = f[i-1] + freq

                                if r < f[i] :
                                        selected_index = i
                                        break

                        selected_allele = alleles[i]
                        if strand == refstrand :
                                if refallele.isupper() :
                                        return selected_allele
                                else :
                                        return selected_allele.lower()

                        else :
                                if refallele.isupper() :
                                        return complement( selected_allele )
                                else :
                                        return complement( selected_allele.lower() )
                else :
                        return refallele

        def polymorph_genome( self, genome ) :
                for linpos, snpinfo_string in self.snp.iteritems() :
                        chr, pos = self.linpos2chromosomepos( linpos )  # @ReservedAssignment
                        direction, alleles, frequencies = self.parse_snp_info( snpinfo_string )
                        refallele = genome.get_nucleotide( chr, pos, direction )

                        #skipping if the original genomic sequence is N...
                        if refallele == 'N' or refallele == 'n' :
                                continue

                        r = random()
                        f = [0.0]*len(frequencies)
                        selected_index = len(frequencies)-1

                        if selected_index < 0 :
                                raise Exception( selected_index )
                        for i, freq in enumerate(frequencies) :
                                if i == 0 :
                                        f[i] = freq
                                else :
                                        f[i] = f[i-1] + freq

                                if r < f[i] :
                                        selected_index = i
                                        break

                        selected_allele = alleles[i]

                        #skipping writing if the reference and polymorph is the same!
                        if selected_allele == refallele.upper() :
                                continue

                        #sometimes the format is not right.
                        if not selected_allele.isalpha() :
                                continue

                        if refallele.isupper() :
                                genome.mutate( chr, pos, selected_allele, direction )
                        else :
                                genome.mutate( chr, pos, selected_allele.lower(), direction )

