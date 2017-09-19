#!/usr/bin/python

import getopt
import sys
import re
import ntpath
import misc_utils as utils

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
	print ("\nUsage:\n\n\t", sys.argv[0], "-L </path/to/loci> -P </path/to/popmap>\n")
	print ("Description:\n")
	print("\tFst_filter.py is a program for estimating Fst for RAD loci and SNPs,\n",\
	"\tand filtering datasets based on Fst and missing data.)


	print("""
Input options:

	-M,--maf	: Input multiple alignment MAF file
	-e,--gff	: GFF file containing annotation information [NOT WORKING YET]
	-L,--loci	: For RAD-data, as the \".loci\" output of pyRAD
	-A,--assembly	: Input whole genome assembly as FASTA [NOT WORKING YET]""")


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
			options, remainder = getopt.getopt(sys.argv[1:], 'M:e:L:A:hc:l:t:b:w:Rm:v:n:Ng:GE:D:p:S:F:s:f:QXo:Pk:Kd:T:', \
			["maf=","gff=","loci=","assembly=",'help',"cov=","len=","thresh=",
			"bait=","win_shift=","mult_reg","min_mult=","var_max=","numN=",
			"callN","numG=","callG","gff_type=","dist_r=","tile_min=",
			"select_r=","filter_r=", "threads=", "blastdb=", "fastadb=",
			"select_b=","filter_b=","quiet","expand","out=",
			"plot_all","mask=","no_mask", "flank_dist=","vsearch=",
			"vthreads=","hacker=", "e_value=", "gapopen=", "gapextend=",
			"word_size=", "megablast", "blastn=", "makedb=", "gap_extend=",
			"word=", "mega", "gap_open=", "blast_db=", "fasta_db=", "--wordsize="])
		except getopt.GetoptError as err:
			print(err)
			display_help("\nExiting because getopt returned non-zero exit status.")
			sys.exit(2)
		#Default values for params
		call_help=0 #boolean

		#Input params
		self.alignment=None
		self.gff=None
		self.loci=None
		self.assembly=None

		#Locus filtering params
		self.cov=1
		self.minlen=None
		self.thresh=0.1
		self.mask=0.1

		#Bait params
		self.blen=80
		self.win_width=None
		self.win_shift=1
		self.mult_reg=0 #boolean
		self.min_mult=None
		self.var_max=0
		self.numN=0
		self.callN=0 #boolean
		self.numG=0
		self.callG=0 #boolean
		self.anchor=None

		#target region options
		self.dist_r=None
		self.flank_dist=500
		self.select_r="rand"
		self.filter_r=0 #bool
		self.filter_t_whole=None
		self.filter_r_objects=[]

		#VSEARCH options - deduplication
		self.vsearch = None
		self.vthreads = 4

		#BLAST options - contaminant removal and filtering by specificity
		self.blastdb=None
		self.fastadb=None
		self.evalue = 0.000001
		self.gapopen = 5
		self.gapextend=2
		self.word_size = None
		self.blast_method = "blastn"
		self.blastn = None
		self.makedb = None

		#Bait selection options
		self.overlap=None
		self.bait_shift=None
		self.select_b="tile"
		self.select_b_num=None
		self.filter_b=0 #bool
		self.filter_b_whole=None
		self.filter_b_objects=[]

		#Running options/ shortcuts
		self.no_mask = 0
		self.stfu = 0

		#Output options
		self.expand = 0
		self.out = None
		self.workdir = None
		self.threads = 1

		self.ploidy=2
		self.db="./mrbait.sqlite"

		#HACKER ONLY OPTIONS
		self._noGraph = 0
		self._noWeightGraph = 0
		self._weightMax = 50000 #maximum size to attempt weighted edge resolution
		self._weightByMin = 0
		self._os = None



		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				display_help("Exiting because help menu was called.")
				sys.exit(0)

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			if opt in ('-M', '--maf'):
				self.alignment = arg
			elif opt in ('-h', '--help'):
				pass
			#Input params
			elif opt in ('-e', '--gff'):
				self.gff = arg
			elif opt in ('-L', '--loci'):
				self.loci = arg
			elif opt in ('-A', '--assembly'):
				self.assembly = arg

			#Locus filtering params
			elif opt in ('-c', '--cov'):
				self.cov = int(arg)
			elif opt in ('-l', '--len'):
				self.minlen = int(arg)
			elif opt in ('-t', '--thresh'):
				self.thresh = float(arg)
			elif opt in ('-k', '--mask'):
				self.mask = float(arg)

			#Bait general params
			elif opt in ('-b', '--bait'):
				self.blen = int(arg)
			elif opt in ('-w', '--win_shift'):
				self.win_shift = int(arg)
			elif opt in ('-R', '--mult_reg'):
				self.mult_reg = 1
			elif opt in ('-m', '--min_mult'):
				self.min_mult = int(arg)
			elif opt in ('-v', '--var_max'):
				self.var_max = int(arg)
			elif opt in ('-n', '--numN'):
				self.numN = int(arg)
			elif opt in ('-N', '--callN'):
				self.callN = 1
			elif opt in ('-g', '--numG'):
				self.numG = int(arg)
			elif opt in ('-G', '--callG'):
				self.callG = 1
			elif opt in ('-E', '--gff_type'):
				self.anchor = arg

			#target region opts
			elif opt in ('-D', '--dist_r'):
				self.dist_r = int(arg)
			elif opt in ('-p', '--tile_min'):
				self.tile_min = int(arg)
				self.tiling = 1
			elif opt in ('-d', '--flank_dist'):
				self.flank_dist = int(arg)
				assert isinstance(self.flank_dist, int), "<--flank_dist> must be an integer"
				assert self.flank_dist >= 0, "<--flank_dist> must be an integer greater than zero!"
			elif opt in ('-S', '--select_r'):
				temp = arg.split('=')
				assert len(temp) == 1, "Invalid specification for <--select_r>: %s"%arg
				self.select_r = (temp[0]).lower()
				chars = (['snp','bad','cons','rand'])
				if self.select_r not in chars:
					raise ValueError("Invalid option \"%r\" for <--select_r>" % self.select_r)
			elif opt in ('-F', '--filter_r'):
				self.filter_r = 1 #turn on region filtering
				#temp = arg.split('/') #parse region filtering options
				self.filter_r_whole = arg
				#for sub in temp:
				subopts = re.split('=|,',arg)
				if subopts[0] in ('snp','mask','gc','len', "pw", "blast_n", "blast_x", "blast_a"): #TODO: Add blast options
					assert len(subopts) == 3, "Incorrect specification of option %r for <--filter_r>" %subopts[0]
					if subopts[0] in ('gc', 'mask'):
						assert 0.0 <= float(subopts[1]) < float(subopts[2]) <= 1.0, "In <--filter_r> suboption \"%s\": Min must be less than max"%subopts[0]
						self.filter_r_objects.append(subArg(subopts[0],float(subopts[1]),float(subopts[2])))
					elif subopts[0] == "len":
						assert subopts[1] < subopts[2], "In <--filter_r> suboption \"%s\": Min must be less than max"%subopts[0]
						self.filter_r_objects.append(subArg(subopts[0],int(subopts[1]),int(subopts[2])))
					elif subopts[0] in ("pw", "blast_x", "blast_i", "blast_a"):
						assert 0.0 <= float(subopts[1]) <= 1.0, "In <--filter_r> suboption \"%s\": Values must be given as proportions! "%subopts[0]
						assert 0.0 <= float(subopts[2]) <= 1.0, "In <--filter_r> suboption \"%s\": Values must be given as proportions! "%subopts[0]
						self.filter_r_objects.append(subArg(subopts[0],float(subopts[1]),float(subopts[2])))
					elif subopts[0] == "snp":
						minS = int(subopts[1])
						maxS = int(subopts[2])
						assert 0 <= minS < maxS, "In <--filter_r> suboption \"%s\": Min must be less than max"%subopts[0]
						assert isinstance(minS, int), "--filter_r suboption snp: Minimum must be an integer"
						assert isinstance(maxS, int), "--filter_r suboption snp: Maximum must be an integer"
						self.filter_r_objects.append(subArg(subopts[0],minS,maxS))
				elif subopts[0] in ('rand','gap','bad'):
					assert len(subopts) == 2, "Incorrect specification of option %r for <--filter_r>" %subopts[0]
					self.filter_r_objects.append(subArg(subopts[0],int(subopts[1])))
				else:
					bad_opts("Invalid option %r for <--filter_r>!" %subopts[0])

			#Bait selection options
			elif opt in ('-s', '--select_b'):
				subopts = re.split('=|,',arg)
				self.select_b = (subopts[0]).lower()
				chars = (['tile', 'center', 'flank'])
				if self.select_b not in chars:
					raise ValueError("Invalid option \"%r\" for <--select_b>" % self.select_b)
				subchars = (['center','flank'])
				if self.select_b in subchars:
					assert len(subopts) == 3, "Incorrect specification of option %r for <--select_b>" %subopts[0]
					self.overlap = int(subopts[2])
					self.select_b_num = int(subopts[1])
					assert self.overlap < self.blen, "Overlap distance cannot be greater than bait length"
				elif self.select_b == "tile":
					self.select_b_num = None
					self.overlap = int(subopts[1])
				#print("select_b is %r" %self.select_b)
				#print("select_b_dist is %r"%self.select_b_dist)
			elif opt in ('-f', '--filter_b'):
				self.filter_b = 1 #turn on region filtering
				#temp = arg.split('/') #parse region filtering options
				self.filter_b_whole = arg
				#for sub in temp:
				subopts = re.split('=|,',arg)
				if subopts[0] in ('min','max','mask','gc'):
					assert len(subopts) == 3, "Incorrect specification of option %r for <--filter_b>" %subopts[0]
					if subopts[0] in ('gc', 'mask'):
						assert subopts[1] < subopts[2], "In <--filter_b> for suboptions \"mask\" and \"gc\": Min must be less than max"
						self.filter_r_objects.append(subArg(subopts[0],float(subopts[1]),float(subopts[2])))
					else:
						self.filter_r_objects.append(subArg(subopts[0],int(subopts[1]),int(subopts[2])))
				elif (subopts[0] is 'rand'):
					assert len(subopts) == 2, "Incorrect specification of option %r for <--filter_b>" %subopts[0]
					self.filter_b_objects.append(subArg(subopts[0],subopts[1]))
				else:
					bad_opts("Invalid option %r for <--filter_b>!" %subopts[0])

			#Running options
			elif opt in ('-Q', '--quiet'):
				self.stfu = 1
			elif opt in ('-K', '--no_mask'):
				self.no_mask = 1

			#vsearch options
			elif opt == ("--vsearch"):
				self.vsearch = str(arg)
			elif opt == ("--vthreads"):
				self.vthreads = int(arg)

			#BLAST options
			elif opt in ("--blastdb", "--blast_db"):
				self.blastdb = arg
			elif opt in ("--fastadb", "--fasta_db"):
				self.fastadb = arg
			elif opt in ("--e_value", "--evalue"):
				self.evalue = float(arg)
			elif opt in ("--gapopen", "--gap_open"):
				self.gapopen = int(arg)
			elif opt in ("--gapextend", "--gap_extend"):
				self.gapextend = int(arg)
			elif opt in ("--word_size", "--word", "--wordsize"):
				self.word_size = int(arg)
			elif opt in ("--megablast", "--mega"):
				self.blast_method = "megablast"
			elif opt == "--blastn":
				self.blastn = arg
			elif opt == "--makedb":
				self.makedb = arg

			#output options
			elif opt in ('-X', '--expand'):
				self.expand = 1
			elif opt in ('-o', '--out'):
				self.out = arg

			#HACKER ONLY OPTIONS
			elif opt in ('--hacker'):
				print(opt, arg)
				subopts = re.split('=|, ',arg)
				main = subopts[0]
				print("Warning: Using \"hacker only\" option <%s>. Be careful."%main)
				if main == "noGraph":
					self._noGraph = 1
				elif main == "win_width":
					assert len(subopts) == 2, "Warning: HACKER option <win_width> must have two arguments separated by \"=\""
					win_width = int(subopts[1])
				elif main == "noWeightGraph":
					self._noWeightGraph = 1
				elif main == "weightByMin":
					self._weightByMin = 1
				elif main == "weightMax":
					assert len(subopts) == 2, "Warning: HACKER option <graphMax> must have two arguments separated by \"=\""
					self._weightMax = int(subopts[1])
				elif main == "os":
					assert len(subopts) == 2, "Warning: HACKER option <os> must have two arguments separated by \"=\""
					os_in = subopts[1]
					if os_in.lower() in ("linux", "ubuntu", "unknown"):
						self._os = "linux"
					elif os_in.lower() in ("mac", "apple", "macos", "darwin", "unix"):
						self_os = "darwin"
					else:
						assert False, "Unrecognized option %s for <--hacker os=>"%os_in
				else:
					assert False, "Unhandled option %r"%main
			else:
				assert False, "Unhandled option %r"%opt

		#DEBUG PRINTS
		for subopt in self.filter_r_objects:
			print("filter_r: Suboption %s has parameters: %s %s" %(subopt.o1,subopt.o2,subopt.o3))

		#for subopt in self.filter_b_objects:
			#print("filter_b: Suboption %s has parameters: %s %s" %(subopt.o1,subopt.o2,subopt.o3))

		#Assertions and conditional changes to params
		if (self.alignment is None) and (self.loci is None) and (self.assembly is None):
			display_help("Input not specified!")
			sys.exit(0)

		assert self.blen > 0, "Bait length cannot be less than or equal to zero!"

		#Set default win_width
		if self.win_width is None:
			self.win_width = self.blen

		#Warn of argument masking
		if self.mult_reg is 0:
			if self.dist_r is not None:
				print("Warning: You have set <--dist_r/-D>, but it is ignored when <--mult_reg/-R> is off. ")

		#Default of dist_r
		if self.dist_r is None:
			self.dist_r = 100

		#Default of overlap
		if self.overlap is None:
			self.overlap = (self.blen // 2)
		self.bait_shift = self.blen - self.overlap

		#if --no_mask, set mask thresh to 1.0
		if self.no_mask:
			self.mask = 1.0

		#Get working dir path and output prefix
		if self.out is None:
			self.out = "mrbait"
			self.workdir = utils.getWorkingDir()
		else:
			self.workdir, self.out = ntpath.split(self.out)
			if self.out == "":
				self.out = "mrbait"
			if self.workdir == "":
				self.workdir = os.getcwd()
		print("Working directory: ", self.workdir)
		print("Prefix is: ", self.out)


		#Assert that win_shift cannot be larger than blen
		if self.minlen is None:
			self.minlen = self.blen
		elif self.blen > self.minlen:
			self.minlen = self.blen

		#set default of min_mult
		if self.min_mult is None:
			self.min_mult = self.minlen

		#If vsearch path not given, try to figure it out
		if self.vsearch is None:
			os_platform = utils.getOS()
			#print("Found OS platform:", os_platform)
			if os_platform == "linux" or os_platform == "unknown":
				print("Automatically detected LINUX or UNKNOWN platform: Using LINUX VSEARCH executable.")
				self.vsearch = utils.getScriptPath() + "/bin/vsearch-2.4.4-linux"
			elif os_platform == "darwin": #mac os
				print("Automatically detected MACOS platform: Using MACOS VSEARCH executable.")
				self.vsearch = utils.getScriptPath() + "/bin/vsearch-2.4.4-macos"

		#BLAST defaults
		if self.word_size is None:
			if self.blast_method == "megablast":
				self.word_size = 28
			else:
				self.word_size = 11

		if self.blastn is None:
			os_platform = utils.getOS()
			#print("Found OS platform:", os_platform)
			if os_platform == "linux" or os_platform == "unknown":
				print("Automatically detected LINUX or UNKNOWN platform: Using LINUX BLASTN executable.")
				self.blastn = utils.getScriptPath() + "/bin/ncbi-blastn-2.6.0-linux"
			elif os_platform == "darwin": #mac os
				print("Automatically detected MACOS platform: Using MACOS BLASTN executable.")
				self.blastn = utils.getScriptPath() + "/bin/ncbi-blastn-2.6.0-macos"

		if self.makedb is None:
			os_platform = utils.getOS()
			#print("Found OS platform:", os_platform)
			if os_platform.lower() in ("linux", "ubuntu", "unknown", None):
				print("Automatically detected LINUX or UNKNOWN platform: Using LINUX BLASTN executable.")
				self.makedb = utils.getScriptPath() + "/bin/ncbi-makeblastdb-2.6.0-linux"
			elif os_platform.lower() in ("darwin", "mac", "macos", "unix", "apple"): #mac os
				print("Automatically detected MACOS platform: Using MACOS MAKEBLASTDB executable.")
				self.makedb = utils.getScriptPath() + "/bin/ncbi-makeblastdb-2.6.0-macos"

		#Default bait design behavior
		if self.select_b == "tile" and self.overlap is None:
			self.overlap = self.blen // 2
