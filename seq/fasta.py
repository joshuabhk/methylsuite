'''
Created on Aug 26, 2014

@author: Bong-Hyun Kim
'''
import cStringIO
from sys import stderr

class FASTA(object):
    '''
    Simple utility class for FASTA file reading and writing.
    '''
    def __init__(self, fp=None, verbose=False):
        '''
        Constructor
        fp: file pointer
        '''
        self.fp = fp
        self.verbose = verbose

        
    def get_record(self):
        fp = self.fp
        
        if self.verbose :
            print >>stderr, 'reading FASTA file:', fp
        
        line = fp.readline()
        
        if self.verbose :
            print >>stderr, 'first line:', line.strip()
            
        while not line.startswith('>') :
            line = fp.readline()
            if not line :
                return None, None
        header = line.strip()
        
        if self.verbose :
            print >>stderr, 'header:', header
            
        line = fp.readline()
        sfp = cStringIO.StringIO()
        while not line.startswith('>'):
            if not line :
                break
            sfp.write(line.strip())
            line = fp.readline()
            
            #if self.verbose :
            #    print >>stderr, 'line:', line.strip()
                
        else :
            fp.seek(-len(line), 1) #rewind if a FASTA heaader appears
            
        return header, sfp.getvalue()
    
        
    def save_as_FASTA(self, fp, name, sequence):
        print >> fp, ">"+name
        print >> fp, sequence
        