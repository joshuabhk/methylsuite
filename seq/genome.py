import cStringIO
from . import reverse_complement, chg_regex, chh_regex, complement
from random import getrandbits
import tempfile
from fasta import FASTA
from sys import stderr


class DuplicateChromosomeError(Exception):
    pass


class NoChromosomeFoundError(Exception):
    pass


class ChromosomeFASTAFromatError(Exception):
    pass


class NotCytosineError(Exception):
    pass


class OutOfBoundChromosomeError(Exception):
    pass


class Chromosome:
    def __init__(self, save_memory=False, tempdir=None, fp=None, verbose=False):
        self.name = None
        self.length = 0
        self.save_memory = save_memory
        self.tempdir = tempdir
        self.fp = None
        self.verbose = verbose

        if self.verbose:
            print >> stderr, "Initializing Chromosome"

        if fp:
            self.read_file(fp)

    def __len__(self):
        return self.length

    def read_file(self, fp):
        '''
        Read the first FASTA record from the content of fp,
        and set the chromosome name and sequence using set_chromosome method.
        '''
        if self.verbose:
            print >> stderr, "reading a FASTA record to set a chromosome"
        fasta = FASTA(fp=fp, verbose=self.verbose)
        chr_name, chr_seq = fasta.get_record()

        if chr_name and chr_seq:
            chr_name = chr_name[1:]
            self.set_chromosome(chr_name, chr_seq)
        elif not chr_name and not chr_seq:
            raise NoChromosomeFoundError(fp.name, chr_name, chr_seq)
        else:
            raise ChromosomeFASTAFromatError(fp.name, chr_name, chr_seq)


    def set_chromosome(self, name, seq):
        if self.save_memory:
            self.fp = tempfile.TemporaryFile(dir=self.tempdir)
        else:
            self.fp = cStringIO.StringIO()

        self.fp.write(seq)
        self.length = len(seq)
        self.fp.seek(0)
        self.name = name

        if self.verbose:
            print >> stderr, self.name, self.length

    def itermethylpositions(self):
        # iteration for call 'C' or 'G' positions.
        fp = self.fp  # reset the reading point
        fp.seek(0)
        i = 0
        s = fp.read(1)
        while s:
            if s == 'C' or s == 'c' or s == 'G' or s == 'g':
                yield i
            s = fp.read(1)
            i += 1

    def itermethylpositions2(self):
        # iteration for call 'C' or 'G' positions.
        # Unlike itermethylpositions,
        # returns chromosome name and position.
        fp = self.fp  # reset the reading point
        fp.seek(0)
        i = 0
        s = fp.read(1)
        while s:
            if s == 'C' or s == 'c' or s == 'G' or s == 'g':
                yield self.name, i
            s = fp.read(1)
            i += 1

    def itermethylposseq(self):
        # iteration for call 'C' or 'G' positions.
        sfp = self.fp
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        while s:
            if s == 'C' or s == 'c' or s == 'G' or s == 'g':
                yield i, s
            s = sfp.read(1)
            i += 1

    def itermethylposseq2(self):
        # iteration for call 'C' or 'G' positions.
        sfp = self.fp
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        while s:
            if s == 'C' or s == 'c' or s == 'G' or s == 'g':
                yield self.name, i, s
            s = sfp.read(1)
            i += 1

    def iterposnt(self):
        sfp = self.fp
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        while s:
            yield i, s
            s = sfp.read(1)
            i += 1

    def iterposnt2(self):
        sfp = self.fp
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        while s:
            yield self.name, i, s
            s = sfp.read(1)
            i += 1

    def peek_3f(self):
        '''Peek 3 NT forward direction including current NT'''
        self.fp.seek(-1, 1)
        nbytes = min(self.length - self.fp.tell(), 3)
        if nbytes > 1:
            n = self.fp.read(nbytes)
            self.fp.seek(-nbytes + 1, 1)  # return to current pos
        elif nbytes == 1:
            n = self.fp.read(nbytes)
        else:
            assert False

        return n

    def peek_3b(self):
        '''Peek backward direction'''
        nbytes = min(3, self.fp.tell())

        if nbytes > 0:
            self.fp.seek(-nbytes, 1)
            n = self.fp.read(nbytes)
        else:
            assert False
        return reverse_complement(n)

    def determine_methyl_context(self, n):
        if n.startswith('C'):
            if n.startswith('CG'):
                return 'CpG'
            elif n.endswith('G'):
                return 'CHG'
            else:
                return 'CHH'
        else:
            raise NotCytosineError(n)

    def itermethylitems(self):
        # iteration to call 'C' or 'G' positions and methylation context
        # CpG, CHG, CHH
        # iteration for call 'C' or 'G' positions.
        sfp = self.fp
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        determine_context = self.determine_methyl_context
        while s:
            if self.verbose > 3:
                print >> stderr, i, s

            if s == 'C' or s == 'c':
                context = self.peek_3f().upper()
                yield i, '+', determine_context(context)

            elif s == 'G' or s == 'g':
                context = self.peek_3b().upper()
                yield i, '-', determine_context(context)

            s = sfp.read(1)
            i += 1

    def itermethylitems2(self):
        # iteration to call 'C' or 'G' positions and methylation context
        # CpG, CHG, CHH
        # iteration for call 'C' or 'G' positions.
        sfp = self.fp
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        determine_context = self.determine_methyl_context
        while s:
            if s == 'C' or s == 'c':
                context = self.peek_3f().upper()
                yield self.name, i, '+', determine_context(context)

            if s == 'G' or s == 'g':
                context = self.peek_3b().upper()
                yield self.name, i, '-', determine_context(context)

            s = sfp.read(1)
            i += 1

    def get_length(self):
        return self.__len__()

    def get_sequence(self, start, length, direction=1):
        # direction can be 0, 1, 2.
        # 0 means random selection
        # 1 means positive strand
        # 2 means reverse complement strand
        if not direction:
            direction = getrandbits(1) + 1

        if start < self.length and start + length < self.length:
            pass
        else:
            raise OutOfBoundChromosomeError(
                self.name, start, length, self.length)

        strfp = self.fp
        strfp.seek(start)
        if direction == 1 or direction == '+':
            return strfp.read(length), '+'
        else:
            return reverse_complement(strfp.read(length)), '-'

    def get_methylation_context(self, start, direction=1):
        # slightly different from the logic of the same method
        # in the Genome class
        # Genome class had inconsistency between get_methylation_context and
        # the itermethylitems method.
        self.fp.seek(start)
        if direction == 1 or direction == '+':
            seq = self.peek_3f().upper()
        else:
            seq = self.peek_3b().upper()

        return self.determine_methyl_context(seq)


    def is_CpG( self, start, direction=None ) :
	# Test if the given site is CpG
	# if the direction is not give, both direction will be tried.
	if direction != None :
		try :
			return 'CpG' == self.get_methylation_context( start, direction )
		except NotCytosineError :
			return False
	else :
		return self.is_CpG( start, 1 ) or self.is_CpG(start, 2 )


    def is_CHG( self, start, direction=None ) :
	# Test if the given site is CHG
	# if the direction is not give, both direction will be tried.
	if direction != None :
		try :
			return 'CHG' == self.get_methylation_context( start, direction )
		except NotCytosinError :
			return False
	else :
		return self.is_CHG( start, 1 ) or self.is_CHG(start, 2 )


    def is_CHH( self, start, direction=None ) :
	# Test if the given site is CHH
	# if the direction is not give, both direction will be tried.
	if direction != None :
		try :
			return 'CHH' == self.get_methylation_context( start, direction )
		except NotCytosinError :
			return False
	else :
		return self.is_CHH( start, 1 ) or self.is_CHH(start, 2 )


    def get_nucleotide(self, start, direction='+'):
        return self.get_sequence(start, length=1, direction=direction)[0]


    def mutate(self, pos, nucleotide, direction='+'):
        if pos < self.length:
            pass
        else:
            raise OutOfBoundChromosomeError(self.name, pos, self.length)

        strfp = self.fp
        strfp.seek(pos)
        if direction == '+' or direction == 1:
            strfp.write(nucleotide)
        else:
            strfp.write(complement(nucleotide))

    def save(self, fp, width=50, transtable=None):
        seqfp = self.fp
        seqfp.seek(0)
        print >> fp, ">" + self.name
        seq = seqfp.read(width)
        if transtable:
            seq = seq.translate(transtable)
        while seq:
            print >> fp, seq
            seq = seqfp.read(width)
            if transtable:
                seq = seq.translate(transtable)


class Genome2:
    '''
    A new class replacing the Genome class.
    It should have the same public method interface as the Genome class 
    with additional methods added to the Genome2 class.
    '''

    def __init__(self, fa=None, save_memory=False, tempdir=None, verbose=False):
        # setting up FASTA file for extracting chromosomes
        self.chromosomes = {}
        self.chromosome_names = []  # ordering information

        self.save_memory = save_memory
        self.tempdir = tempdir
        self.verbose = verbose
        if fa:
            self.add_chromosomes(fa)

    def __contains__(self, chr_name):
        return chr_name in self.chromosomes

    def __getitem__(self, chr_name):
        return self.chromosomes[chr_name]

    def add_chromosome(self, chromosome):
        if chromosome.name in self.chromosomes:
            raise DuplicateChromosomeError(chromosome.name)
        self.chromosomes[chromosome.name] = chromosome
        self.chromosome_names.append(chromosome.name)

    def add_chromosomes(self, fa):
        if not fa:
            return

        fp = open(fa)
        try:
            chromosome = Chromosome(save_memory=self.save_memory,
                                    tempdir=self.tempdir,
                                    fp=fp, verbose=self.verbose)
            while chromosome:
                self.add_chromosome(chromosome)
                chromosome = Chromosome(save_memory=self.save_memory,
                                        tempdir=self.tempdir,
                                        fp=fp, verbose=self.verbose)
        except NoChromosomeFoundError:
            if self.chromosomes:
                pass
            else:
                raise

    def itermethylpositions(self, chromosome_name=None):
        # iteration for call 'C' or 'G' positions.
        if chromosome_name:
            for i in self.chromosomes[chromosome_name].itermethylpositions():
                yield i
        else:
            for chr_name in self.chromosome_names:
                for i in self.chromosomes[chr_name].itermethylpositions2():
                    yield i

    def itermethylposseq(self, chromosome_name=None):
        # iteration for call 'C' or 'G' positions.
        if chromosome_name:
            for i in self.chromosomes[chromosome_name].itermethylposseq():
                yield i
        else:
            for chr_name in self.chromosome_names:
                for i in self.chromosomes[chr_name].itermethylposseq2():
                    yield i

    def iterposnt(self, chromosome_name=None):
        # iteration for call 'C' or 'G' positions.
        if chromosome_name:
            for i in self.chromosomes[chromosome_name].iterposnt():
                yield i
        else:
            for chr_name in self.chromosome_names:
                for i in self.chromosomes[chr_name].iterposnt2():
                    yield i

    def itermethylitems(self, chromosome_name=None):
        if chromosome_name:
            for i in self.chromosomes[chromosome_name].itermethylitems():
                yield i
        else:
            for chr_name in self.chromosome_names:
                for i in self.chromosomes[chr_name].itermethylitems2():
                    yield i

    def get_length(self, chr_name):
        return len(self.chromosomes[chr_name])

    def get_sequence(self, chromosome_name, start, length, direction=1):
        return self.chromosomes[chromosome_name] \
            .get_sequence(start, length, direction)

    def get_methylation_context(self, chromosome_name, start, direction=1):
        return self.chromosomes[chromosome_name] \
            .get_methylation_context(start, direction)

    def get_nucleotide(self, chromosome_name, start, direction='+'):
        return self.chromosomes[chromosome_name] \
            .get_sequence(start, length=1, direction=direction)

    def mutate(self, chromosome_name, pos, nucleotide, direction='+'):
        self.chromosomes[chromosome_name].mutate(pos, nucleotide,
                                                 direction=direction)

    def save(self, fp, width=50, transtable=None):
        for seqid in sorted(self.chromosome_names):
            self.chromosomes[seqid].save(
                fp, width=width, transtable=transtable)

    def save_chromosome(self, seqid, fp, width=50, transtable=None):
        chromosome = self.chromosomes[seqid]
        chromosome.save(fp, width=width, transtable=transtable)


class Genome:
    def __init__(self, fa):
        # setting up FASTA file for extracting
        # chromosomes
        self.chromosomes = {}
        self.lengths = {}
        if fa:
            fp = open(fa)
            strfp = cStringIO.StringIO()
            for l in fp:
                if l.startswith('>'):
                    # if strfp contains chromosomes
                    if strfp.tell():
                        self.chromosomes[seqid] = strfp  # @UndefinedVariable
                        self.lengths[seqid] = strfp.tell()  # @UndefinedVariable
                        strfp = cStringIO.StringIO()
                    seqid = l[1:].strip()
                else:
                    strfp.write(l.strip())
            else:
                if strfp.tell():
                    self.chromosomes[seqid] = strfp
                    self.lengths[seqid] = strfp.tell()

    def __contains__(self, chromosome):
        return chromosome in self.chromosomes

    def itermethylpositions(self, chromosome):
        # iteration for call 'C' or 'G' positions.
        sfp = self.chromosomes[chromosome]
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        while s:
            if s == 'C' or s == 'c' or s == 'G' or s == 'g':
                yield i
            s = sfp.read(1)
            i += 1

    def itermethylposseq(self, chromosome):
        # iteration for call 'C' or 'G' positions.
        sfp = self.chromosomes[chromosome]
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        while s:
            if s == 'C' or s == 'c' or s == 'G' or s == 'g':
                yield i, s
            s = sfp.read(1)
            i += 1

    def iterposnt(self, chromosome):
        sfp = self.chromosomes[chromosome]
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        while s:
            yield i, s
            s = sfp.read(1)
            i += 1

    def peek_3f(self, fp):
        '''Peek forward direction'''
        fp.seek(-1, 1)
        n = fp.read(3)
        fp.seek(-2, 1)  # return to current pos
        return n

    def peek_3b(self, fp):
        '''Peek backward direction'''
        fp.seek(-3, 1)
        n = fp.read(3)
        return reverse_complement(n)

    def determine_methyl_context(self, n):
        if n.startswith('C'):
            if n.startswith('CG'):
                return 'CpG'
            elif n.endswith('G'):
                return 'CHG'
            else:
                return 'CHH'
        else:
            raise Exception(
                'Methylation context determination should start from C', n)

    def itermethylitems(self, chrom):
        # iteration to call 'C' or 'G' positions and methylation context
        # CpG, CHG, CHH
        # iteration for call 'C' or 'G' positions.
        sfp = self.chromosomes[chrom]
        sfp.seek(0)  # reset the reading point
        i = 0
        s = sfp.read(1)
        determine_context = self.determine_methyl_context
        while s:
            if s == 'C' or s == 'c':
                context = self.peek_3f(sfp).upper()
                yield i, '+', determine_context(context)

            if s == 'G' or s == 'g':
                context = self.peek_3b(sfp).upper()
                yield i, '-', determine_context(context)

            s = sfp.read(1)
            i += 1

    def get_length(self, chrom):
        return self.lengths[chrom]

    def get_sequence(self, chrom, start, length, direction=1):
        # direction can be 0, 1, 2.
        # 0 means random selection
        # 1 means positive strand
        # 2 means negative strand
        if not direction:
            direction = getrandbits(1) + 1

        if start < self.lengths[chrom] and start + length < self.lengths[chrom]:
            pass
        else:
            raise Exception(
                "OutOfBoundReading", start, length, self.lengths[chrom])

        strfp = self.chromosomes[chrom]
        strfp.seek(start)
        if direction == 1 or direction == '+':
            return strfp.read(length), '+'
        else:
            return reverse_complement(strfp.read(length)), '-'

    def get_methylation_context(self, chrom, start, direction=1):

        if direction == 1 or direction == '+':
            seq = self.get_sequence(chrom, start, 3, direction)[0]
        else:
            seq = self.get_sequence(chrom, start - 2, 3, direction)[0]

        if not seq.isupper():
            seq = seq.upper()

        if seq.startswith('CG'):
            return 'CpG'
        elif chg_regex.match(seq):
            return 'CHG'
        elif chh_regex.match(seq):
            return 'CHH'
        else:
            raise Exception(
                "Methylation context cannot be determined", chrom, start, direction, seq)

    def get_nucleotide(self, chrom, start, direction='+'):
        strfp = self.chromosomes[chrom]
        strfp.seek(start)

        if start < self.lengths[chrom]:
            pass
        else:
            raise Exception("OutOfBoundReading", start, self.lengths[chrom])

        if direction == '+' or direction == 1:
            return strfp.read(1)
        else:
            return complement(strfp.read(1))

    def mutate(self, chrom, pos, nucleotide, direction='+'):
        # assumes + direction
        if pos < self.lengths[chrom]:
            pass
        else:
            raise Exception("OutOfBound", chrom, pos, self.lengths[chrom])

        strfp = self.chromosomes[chrom]
        strfp.seek(pos)
        if direction == '+' or direction == 1:
            strfp.write(nucleotide)
        else:
            strfp.write(complement(nucleotide))

    def save(self, fp, width=50, transtable=None):
        for seqid in sorted(self.chromosomes):
            self.save_chromosome(seqid, fp, width=width, transtable=transtable)

    def save_chromosome(self, seqid, fp, width=50, transtable=None):
        seqfp = self.chromosomes[seqid]
        seqfp.seek(0)
        print >> fp, ">" + seqid
        seq = seqfp.read(width)
        if transtable:
            seq = seq.translate(transtable)
        while seq:
            print >> fp, seq
            seq = seqfp.read(width)
            if transtable:
                seq = seq.translate(transtable)
