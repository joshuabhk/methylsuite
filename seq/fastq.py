import bz2
from statlib.stats import mean, stdev
import cStringIO


def quality_string2score(quality_string, base=33):
    # assuming Sanger sequencing quality score #base 33
    return [ord(a) - base for a in quality_string]


class FASTQ:

    def __init__(self, fn=None, fp=None):
        self.fastq = ''
        self.quality_strings = []
        self.fp = None
        if fn:
            self.fastq = fn
            if fn.endswith('.bz2'):
                self.fp = bz2.BZ2File(fn)
            else:
                self.fp = open(self.fastq)
        elif fp:
            self.fp = fp

    def reset(self):
        self.fp.seek(0)

    def read_all_quality_strings_by_position(self):
        currentfp = self.fp.tell()  # mark current file pointer
        self.fp.seek(0)  # rewind just in case

        rec = self.get_record()
        # use the length of first record
        # assuming that the length of sequences are all same
        length = len(rec[-1])
        strfps = []
        for i in xrange(length):
            strfps.append(cStringIO.StringIO())

        while rec:
            [strfps[i].write(rec[-1][i]) for i in xrange(length)]
            rec = self.get_record()
        else:
            self.quality_strings = [strfp.getvalue() for strfp in strfps]

        self.fp.seek(currentfp)  # rewind after usage

    def set_quality_mean_and_variation_per_position(self):
        self.read_all_quality_strings_by_position()
        self.means = [mean(self.quality_string2score(s))
                      for s in self.quality_strings]
        self.stdevs = [stdev(self.quality_string2score(s))
                       for s in self.quality_strings]

        # relase the memory for quality_strings
        self.quality_strings = None

    def force_sanger_range(self, n):
        # return value >33 and <73.
        if n < 33:
            return 33
        elif n > 73:
            return 73
        else:
            return n

    def build_high_quality_string(self, length=82, base=33):
        scores = [chr(34 + base)] * length
        return ''.join(scores)

    def build_low_quality_string(self, length=82, base=33):
        scores = [chr(2 + base)] * length
        return ''.join(scores)

    def get_record(self, autorewind=False):
        if self.fp:
            header1 = self.fp.readline().strip()
            seq = self.fp.readline().strip()
            header2 = self.fp.readline().strip()
            qual = self.fp.readline().strip()

            if autorewind:
                if header1 == seq == header2 == qual == '':
                    self.fp.seek(0)
                    header1 = self.fp.readline().strip()
                    seq = self.fp.readline().strip()
                    header2 = self.fp.readline().strip()
                    qual = self.fp.readline().strip()
            else:
                if header1 == seq == header2 == qual == '':
                    return None

            if not header1.startswith('@'):
                raise Exception(self.fn)

            if not header2.startswith('+'):
                raise Exception(self.fn)

            return (header1, seq, header2, qual)

    def determine_trim(self, qual, cutoff, base=33):
        scores = self.quality_string2score(qual, base=base)
        for i, s in enumerate(scores):
            if s <= cutoff:
                return i
        return len(scores)

    def write_record(self, sfp, header, seq, header2, qual, 
                     trimcutoff=None, base=33, minlen=0, transtable=None):
        # This function will truncate the read at the position where the
        # trimcutoff reached!

        trim = self.determine_trim(qual, trimcutoff, base=base)
        if trim == 0:
            return 0
        if trim < minlen:
            return 0
        else:
            seq = seq[:trim]
            qual = qual[:trim]

        print >>sfp, header

        if transtable:
            print >>sfp, seq.translate(transtable)
        else:
            print >>sfp, seq

        print >>sfp, header2
        print >>sfp, qual
        return 1

    def iterreadids(self):
        r = self.get_record()
        while r:
            yield r[0].replace('@', '')
            r = self.get_record()

    def iterrecords(self):
        r = self.get_record()
        while r:
            yield r
            r = self.get_record()

    def save(self, fp, transtable=None):
        '''
        For faster performance, it is implemented not using
        iterreocrds/write functions in the class.
        For speed, I streamlined read/write with translation function.
        '''
        self.reset()
        for i, line in enumerate(self.fp):
            if i % 4 == 1:
                if transtable:
                    line = line.translate(transtable)
                fp.write(line)
            else:
                fp.write(line)

    def get_next_quality_score_line(self, autorewind=True):
        h1, s, h2, q = self.get_record(autorewind=autorewind)  # @UnusedVariable
        return q

    def quality_string2score(self, quality_string, base=33):
        # assuming Sanger sequencing quality score #base 33
        return [ord(a) - base for a in quality_string]

    def score2quality_string(self, scores, base=33):
        return ''.join([chr(i + base) for i in scores])

    def quality_score2error_rate(self, Q):
        return pow(10, -Q / 10.0)
