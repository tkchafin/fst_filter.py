#!/usr/bin/python

import getopt
import sys

def string_containsAny(str, set):
	for c in set:
		if c in str: return 0;
	return 1;

def bad_opts(message=""):
	print(message)
	print("Invalid options: Exiting program. Please see manual or use <--help>\n")
	sys.exit(1)

def display_help(message=None):
	if message is not None:
		print (message)
	print ("\nFst_filter.py\n")
	print ("Contact:\n\n\tTyler K. Chafin\n\tUniversity of Arkansas\n\ttkchafin@uark.edu\n")
	print ("\nUsage:\n\n\t", sys.argv[0], "-i </path/to/phylip> -p </path/to/popmap>\n")
	print ("Description:\n")
	print("\tFst_filter.py is a program for estimating Fst for RAD loci and SNPs,\n",\
	"\tand filtering datasets based on Fst and missing data.")


	print("""
Input options:

	-i,--input	: Input SNP alignment (PHYLIP-formatted)
	-p,--popmap	: Tab-delimited popmap file

""")


#Sub-object for holding filtering options
class subArg():
	def __init__(self, o1, o2, o3=None):
		self.o1 = o1
		self.o2 = o2
		self.o3 = o3

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'i:p:h', \
			["help","input=","popmap="])
		except getopt.GetoptError as err:
			print(err)
			display_help("\nExiting because getopt returned non-zero exit status.")
			sys.exit(2)

		#Input params
		self.input=None
		self.popmap=None

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				display_help("Exiting because help menu was called.")
				sys.exit(0)

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			if opt in ('-i', '--input'):
				self.input = arg
			elif opt in ('-p', '--popmap'):
				self.popmap = arg
			elif opt in ('-h', '--help'):
				pass
			else:
				assert False, "Unhandled option %r"%opt

		if self.input is None:
			display_help("SNP file not defined!")
			sys.exit(1)

		if self.popmap is None:
			display_help("POPMAP file not defined!")
			sys.exit(1)
