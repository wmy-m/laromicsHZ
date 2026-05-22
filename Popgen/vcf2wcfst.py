#! /usr/bin/env python

# Author: (c) David Marques, 07.11.2024
# Title: vcf2wcfst.py
# Written in Python 3.11.5

# What it does: This script computes Weir & Cockerham's FST 
# (Weir & Cockerham 1984 Evolution) from a VCF file with genotypes
# and two population files (list of individuals, one per line),
# and outputs not only FST, but also the Numerator (N, a in Weir & Cockerham)
# and Denominator (D, a+b+c in Weir & Cockerham), allowing to recalculate 
# weighted mean FST across multiple output files by averaging N's and D's
# across chunks and then building the ratio (see Bhatia et al. 2013)

# Tested against implementation in vcftools 0.1.16 -> same values output

# Import python modules
from sys import *
import argparse, re
import gzip

# Fixes a broken pipe issue if STDOUT is piped into head
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

# Read input arguments using the argparse module
parser = argparse.ArgumentParser(description='Computes Weir & Cockerhams FST from a VCF file and two population files (lists of individuals, one per line)')
parser.add_argument('-i', '--input', dest='i', help="Input file in VCF format or enter '-' for STDIN [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="Prefix for output file or enter '-' for STDOUT [required]", required=True)
parser.add_argument('-p', '--popfile1', dest='p', help="Population file 1 [required], text file with individuals listed one per line", required=True)
parser.add_argument('-q', '--popfile2', dest='q', help="Population file 2 [required], text file with individuals listed one per line", required=True)
parser.add_argument('-m', '--minGT', dest='m', help="Sets the minimum number of genotypes per population for an FST to be computed. Optional, default and lowest possible value is 1.", required=False, default=1)

args = parser.parse_args()

# Read population files and store information in vector
inputP=open(args.p,'r')
indpop1=[]
for Line in inputP:
	indpop1+=[Line.strip("\n")]
inputP.close()
inputQ=open(args.q,'r')
indpop2=[]
for Line in inputQ:
	indpop2+=[Line.strip("\n")]
inputQ.close()

# Weir and Cockerham's FST Function
# Input: gt1 and gt2 are arrays of 0, 1, 2 genotype format, thus only
# containing non-missing genotypes at biallelic SNPS
def wcFST(i1,i2): # Weir & Cockerhams FST function
	# Define the number of populations
	r=2
	# Compute n1, n2 number of diploid individuals (sample size) per population
	n1=len(i1)
	n2=len(i2)
	# Combute n_bar, the average sample size
	nbar=(n1/r)+(n2/r)
	if nbar > 1: # Otherwise leads to division by zero
		# Compute nc the squared coefficient of variaten of sample sizes
		nc=(r*nbar-(((n1**2)/(r*nbar))+((n2**2)/(r*nbar))))/(r-1)
		# Compute p1, p2 allele frequency in populations 1 and 2
		p1=(int(sum(i1))/n1)/2
		p2=(int(sum(i2))/n2)/2
		# Compute p_bar, the average sample frequency of allelel A
		pbar=(n1*p1/(r*nbar))+(n2*p2/(r*nbar))
		# Compute s_squared, the sample variance of allele frequencies over populations
		s2=((n1*((p1-pbar)**2))/((r-1)*nbar))+((n2*((p2-pbar)**2))/((r-1)*nbar))
		# Compute h1 and h2, heterozygote frequencies in populations 1 and 2
		h1=sum([x==1 for x in i1])/n1
		h2=sum([x==1 for x in i2])/n2
		# Compute h_bar, the average heterozygote frequency for allele A
		hbar=((n1*h1)/(r*nbar))+((n2*h2)/(r*nbar))
		# Compute components a, b and c
		a=nbar/nc*(s2-(1/(nbar-1))*((pbar*(1-pbar))-((r-1)/r*s2)-(1/4*hbar)))
		b=(nbar/(nbar-1))*((pbar*(1-pbar))-((r-1)/r*s2)-((2*nbar-1)/(4*nbar)*hbar))
		c=1/2*hbar
		# Compute Numerator and Denominator
		N=a
		D=(a+b+c)
		fst=N/D
		return str(N)+"\t"+str(D)+"\t"+str(fst)
	else:
		return "-nan\t-nan\t-nan"

# Parse -i option, normal input VCF file, gzipped input VCF file or STDIN
if args.i == "-":
	inputF=stdin
else:
	# Decide whether to open gzipped file or not
	if args.i.endswith(".gz"):
		inputF=gzip.open(args.i,'rt')
	else:
		inputF=open(args.i,'r')

# Decide whether to write into STDOUT or file
if args.o =="-":
	socket=stdout
else:
	outputF=open(args.o+".wcfst.txt", 'w')
	socket=outputF

# Working across VCF file
for Line in inputF:
	# HEADER SECTION of VCF file: parse header information on individuals / populations
	if re.match('^#',Line): 
		if re.match('^##',Line) is None: # header with individuals / popinfo for parsing / changing
			header=Line.strip("\n").split("\t")
			indid=header[9:len(header)]  # header now contains all the individual IDs
			# Find population indices of the two populations
			idxpop1=[indid.index(i) for i in indpop1]
			idxpop2=[indid.index(i) for i in indpop2]
			# Write header
			socket.write("Chr\tPos\tN\tD\tFST\n")
	# DATA SECTION of VCF file: fill SFS entries
	else:
		columns=Line.strip("\n").split("\t")
		# Check if site is biallelic SNP (and not indel etc.)
		if (len(columns[3])==1) and (len(columns[4])==1):
			genotypecolumns=columns[9:len(columns)]
			tmp=[x.split(":") for x in genotypecolumns]
			genotypes=[x[0] for x in tmp]
			gt1=[genotypes[i] for i in idxpop1]
			gt2=[genotypes[i] for i in idxpop2]
			# Remove missing genotypes
			gt1=list(filter(('./.').__ne__, gt1))
			gt2=list(filter(('./.').__ne__, gt2))
			gt1=list(filter(('.').__ne__, gt1))
			gt2=list(filter(('.').__ne__, gt2))
			# Convert genotypes to 0, 1, 2 format
			i1=[int(y)+int(z) for [y,z] in [x.replace("|","/",1).split("/") for x in gt1]]
			i2=[int(y)+int(z) for [y,z] in [x.replace("|","/",1).split("/") for x in gt2]]
			# Check if the site is polymorphic
			if (len(set(i1+i2))>1) or (set(i1+i2)=={1}):
				# Check that at least args.m genotypes per population is present
				if (len(i1)>int(args.m)) and (len(i2)>int(args.m)):
					# Write to output
					socket.write(columns[0]+"\t"+columns[1]+"\t"+wcFST(i1,i2)+"\n")

# Close input file if applicable
if args.i != "-":
	inputF.close()

# Close output file if applicable
if args.o !="-":
	outputF.close()
