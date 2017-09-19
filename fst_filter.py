#!/usr/bin/python

import os
import sys
import Bio
from Bio import AlignIO
import aln_file_tools as aft
import alignment_tools as a
from fst_filter_menu import parseArgs

"""IMPORTANT NOTE:
Pausing development for now because I found an easier way to get what I want.
May resume working on this later. """



#Need:
	#Vector of population sizes
	#Mean population size of each



#For each Locus:
	#For each allele:
		#Calculate heterozygosity per pop
		#Calculate frequency per pop
		#Calculate average het of allele
		#Calculate average freq of allele
		#Calculate sample variance of allele
		#Calculate average pop size and variance

params = parseArgs()

aln = aft.read_phylip(params.input)

for col in a.aln_column_iterator(aln):
	print(col)
