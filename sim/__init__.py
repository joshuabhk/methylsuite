from random import randrange, getrandbits, gauss, random, choice

def random_nucleotide( nuc, nucleotide_pool = ['A','C','G','T'] ) :
        if nuc.islower() :
                nuc = nuc.upper()
                #make sure newnuc is different from nuc
                newnuc = choice( nucleotide_pool )
                while newnuc == nuc :
                        newnuc = choice( nucleotide_pool )
                return newnuc.lower()
        else :
                newnuc = choice( nucleotide_pool )
                while newnuc == nuc :
                        newnuc = choice( nucleotide_pool )
                return newnuc


def choice_with_N( quality_score, nucleotide_pool=['A','C','G','T'] ) :
        #select N or other randomly selected nulceotides.
        #print "choice_with_N",

        if quality_score <= 4 :
                if random() < 0.00185 : return 'N'
                else : return choice(nucleotide_pool)
        elif quality_score == 5 :
                if random() < 0.577 : return 'N'
                else : return choice(nucleotide_pool)
        elif quality_score == 6 :
                if random() < 0.0863 : return 'N'
                else : return choice(nucleotide_pool)
        elif quality_score == 7 :
                if random() < 0.00508 : return 'N'
                else : return choice(nucleotide_pool)
        elif quality_score == 8 :
                if random() < 0.0000163 : return 'N'
                else : return choice(nucleotide_pool)
        else :
                choice( nucleotide_pool )

def random_nucleotide2( nuc, nucleotide_pool = ['A','C','G','T'] ) :
        #introdue N if quality score is less than 9.
        #the N introduction rate is determined by my little experiment.
        if nuc.islower() :
                #make sure newnuc is different from nuc
                newnuc = choice( nucleotide_pool )
                return newnuc.lower()
        else :
                newnuc = choice( nucleotide_pool )
                return newnuc


