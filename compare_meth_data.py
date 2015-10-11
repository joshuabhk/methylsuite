#!/usr/bin/env python
import sys



class EmptyLine( Exception ) :
	pass

class MethLine :
	def __init__(self, fp=None) :
		self.fp = fp
		self.empty = True
		self.chr = None
		self.coord = None
		self.mc = None
		self.nceff = None
		self.nc = None
		self.genotype = None
		self.seq = None
		
		self.count = 0

	def parse_line( self ) :
		l = self.fp.readline()
		self.l = l.strip()
		if not l :
			self.empty = True
			return

		l = l.split()
		if not l :
			self.empty = True 
			return

		self.count += 1
		try :
			self.empty = False
			self.chr = l[0]
			self.coord = int(l[1])
			self.mc = int(l[3])
			self.nceff = float(l[4] )
			self.nc = float(l[5])
			self.genotype = l[6]
			self.seq = l[7]

			total = self.mc + self.nc
			if total :
				self.rawmethyl = self.mc/total
			else :
				self.rawmethyl = None
			self.total = total

			total = self.mc + self.nceff
			if total :
				self.effmethyl = self.mc/total
			else :
				self.effmethyl = None
			self.efftotal = total

		except IndexError :
			self.empty = True
		
	
	def __str__( self ) :
		return self.l

	def location_diff( self, m ) :
		if self.empty or m.empty :
			raise EmptyLine( self, m )

		if self.chr == m.chr :
			return self.coord - m.coord
		elif self.chr < m.chr :
			return -1
		else :
			return 1

	def sig_diff( self, m, mr_cutoff=0.25, cov_cutoff=10, methylrate_type='raw' ) :
		diff = self.methyldiff( m, methylrate_type=methylrate_type )
		if diff == None :
			return False

		if abs(diff) > mr_cutoff and self.total > cov_cutoff and m.total>cov_cutoff :
			return True
		else :
			return False
		

	def methyldiff( self, m, methylrate_type='raw' ) :
		if methylrate_type == 'raw' :
			mr1 = self.rawmethyl
			t1 = self.total 
			
			mr2 = m.rawmethyl
			t2 = m.total

		elif methylrate_type == 'eff' :
			mr1 = self.effmethyl
			t1 = self.efftotal
		
			mr2 = m.effmethyl
			t2 = m.efftotal

		if mr1 != None and mr2 != None :
			return mr1 - mr2
		else :
			return None

def main() :
	fn1 = sys.argv[1]
	fn2 = sys.argv[2]

	fp1 = open( fn1 )
	fp2 = open( fn2 )
	
	m1 = MethLine(fp1)
	m2 = MethLine(fp2)

	m1.parse_line()
	m2.parse_line()

	common_lines = 0
	sig_diff_lines = 0
	while not m1.empty :
		location_diff = m1.location_diff(m2) 

		#print >>sys.stderr, m1
		#print >>sys.stderr, m2

		if location_diff == 0 :
			if m1.sig_diff(m2) :
				print m1 + '\t' + m2
				sig_diff_lines += 1
			common_lines += 1

		elif location_diff > 0 :
			m2.parse_line()
				
		elif location_diff < 0 :
			m1.parse_line()
		

	print >>sys.stderr, "T1:", m1.count, 'T2:',m2.count, 'Common:', common_lines, 'Significant:', sig_diff_lines



if __name__ == '__main__' :
	main()
