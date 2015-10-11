from __init__ import NotYetImplemented
from util import Interval

class ReadInterval(Interval) :
    '''
    class for modeling raw read interval and
    query interval to the interval in the chromosome space
    
    Raw Read space <--> Query mapped space <--> Chromosome space
    '''
    def __init__(self, start=None, end=None, interval=None, aligned_read=None):
        self.start = start #interval supposed to be in Raw Read space
        self.end = end #interval supposed to be defined in the Raw Read Space
        if interval :
            self.start = interval.start
            self.end = interval.end
            
        self.aread = aligned_read # the AlignedRead class instance from pysam
        raise NotYetImplemented()
    
    def __contains__(self, ri):
        if isinstance(ri, int) :
            if self.start <= ri < self.end :
                return True
            else :
                return False
        elif isinstance(ri, Interval) :
            return Interval.__contains__(ri)
        else :
            raise NotYetImplemented()
        
    def ri2qi(self, ri):
        raise NotYetImplemented()
            