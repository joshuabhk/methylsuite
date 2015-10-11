
from string import maketrans
import re

#fast calculation of reverse complement
ComplementTable = maketrans( 'AaGgCcTtUuRrSsYyWwKkMmBbVvDdHhNn', 
							 'TtCcGgAaAaYySsRrWwMmKkVvBbHhDdNn' )
def reverse_complement( seq ) :
	return seq.translate( ComplementTable )[::-1]

def complement( seq ):
	return seq.translate( ComplementTable )

#regular expression for cytosine context
chg_regex = re.compile( 'C[^G]G')
chh_regex = re.compile('C[^G][^G]')

#C to T translation
CtoTtable = maketrans( 'cC', 'tT' )

#G to A translation
GtoAtable = maketrans( 'gG', 'aA' )


