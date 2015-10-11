import sys

#A useful function for parsing reference mapping position information
#from the fastq header line
def parse_header( header ) :
        '''
        Parses FASTQ Header when the fields are concaternated by underscores.
        '''
        fields = header.split('_')
        if len(fields) == 6 :
                junk, id, chr, pos, length, direction = fields  # @ReservedAssignment @UnusedVariable
                return id, chr, pos, length, direction
        elif len(fields) > 6 :
                return fields[1], '_'.join( fields[2:-3] ), fields[-3], fields[-2], fields[-1]
        else :
                print >>sys.stderr, header
                assert 0


def find_FASTA_record( fp ):
    '''
    Skip lines until the start of the next FASTA record.
    returns True if the next FASTA record is found,
    returns False if EOF reached before the next FASTA record.
    '''
    pass

def read_FASTA_record( fp ):
    '''
    '''
    pass



class RangeInclusionError(Exception):
    pass


def overlap( range1, range2 ):
    '''
    Check the overlapping between the two ranges.
    Throws an exception when a range is belong to the other.
    '''
    
    if range1 > range2 :
        range1, range2 = range2, range1
        
    if range1[1] > range2[0] :
        if range1[1] > range2[1] :
            raise RangeInclusionError( range1, range2 )
        return True
    else : 
        return False

def split_middle( range1, range2 ):
    '''
    returns two new ranges that split half of the overlapping region
    and makes them exclusively covered.
    
    Note that this method does not check for the validity of the range1 and range2's 
    overlap.
    '''
    switch = False
    if range1 > range2 :
        switch = True
        range1, range2 = range2, range1
    
    r1 = list(range1)
    r2 = list(range2)
    
    r1[1] = (r1[1]+r2[0])/2
    r2[0] = r1[1]+1
    
    if switch :
        return r2, r1
    else :
        return r1, r2
    
class Interval:
    '''
    utility class defining interval and interval operations
    '''
    def __init__(self, start=None, end=None):
        self.start = start
        self.end = end
        
    def __comp__(self, interval):
        startdiff = interval.start - self.start
        if not startdiff :
            return interval.end - self.end
        return startdiff
        
    def is_included(self, interval):
        '''
        test if the given interval is included to self or
        if self is included to the interval
        '''
        if self < interval :
            if self.end > interval.end :
                return True
        else :
            if interval.end > self.end :
                return True
            
    def overlap(self, interval):
        if interval < self :
            if interval.end > self.start :
                return True
        else :
            if self.end > interval.start :
                return True
        return False

    def split_middle(self, interval):
        r1 = Interval( self.start, self.end )
        r2 = Interval( interval.start, interval.end)
        
        if self < interval :
            r1.end = (r1.end+r2.start)/2
            r2.start = r1.end + 1
        else :
            r2.end = (r2.end + r1.start)/2
            r1.start = r2.end + 1
        
        return r1, r2

    
    def __contains__(self, ri):
        pass
    
    
            