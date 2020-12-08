#!/usr/bin/env python

"""
By Veronica Mixao 
@BSC, Barcelona
vmixao@gmail.com

Objective: Detect crossing-over events in heterozygous regions using HapCut2
Last change: 21/11/2019
"""

import argparse
import os
import numpy

bedtools = "/usr/bin/bedtools"
hapcut = "~/HapCut/HapCUT2/hapcut2/build/HAPCUT2"
extract_hairs = "~/HapCut/HapCUT2/hapcut2/build/extractHAIRS"

def cross_over(bamfile,vcffile,tag,min_var,max_switch):
	print "**************************************"
	print "RUNNING PIPELINE TO DETECT CROSSING-OVER ON " + tag
	print "\n"
	
	#run extract_hairs
	print "Obtaining fragments file for " + tag
	out1 = tag + ".fragments"
	os.system(extract_hairs + " --bam " + bamfile + " --VCF " + vcffile + " --PEonly > " + out1)
	
	#run hapcut2
	print "Running HapCUT2 for " + tag
	out2 = tag + ".hapcut2"
	os.system(hapcut + " --fragments " + out1 + " --vcf " + vcffile + " --output " + out2)
	
	#check hapcut results
	print "Obtaining blocks whith two different genotypes for " + tag
	hapcut_file = open(out2, "r")
	hap = hapcut_file.readlines()
	
	out = open(tag + ".m_" + str(min_var) + ".s_" + str(max_switch) + ".cross_over.bed", "w+")
	
	print "Blocks of interest for " + tag
	
	d = {}
	af = []
	
	for line in hap:
		if "********" not in line:
			if "BLOCK:" in line: #a new block is starting
				l = line.split() 
				start = l[2]
				n_var = l[4]
				phased = l[6]
				length = l[8]
				
				#this is not the first block
				if "01" in d.keys() or "10" in d.keys():
					if "01" in d.keys() and "10" in d.keys():
						d[GT].append(count) #I need to append the last count to the last GT
						if len(d["01"]) > 0 and len(d["10"]) > 0: #both haplotypes were observed in the previous block
							if numpy.amax(d["01"]) >= int(min_var) and numpy.amax(d["10"]) >= int(min_var): #both genotypes present blocks with more than 5 variants phased
								if (switches/2)/float(length) < float(max_switch):
									if numpy.mean(af) < 0.6 and numpy.mean(af) > 0.4:
										print chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + block_ID + "\t" + tag + "\t" + str(numpy.amax(d["01"])) + "\t" + str(numpy.amax(d["10"])) #report the results of the previous block
										print >>out, chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + block_ID + "\t" + tag + "\t" + str(numpy.amax(d["01"])) + "\t" + str(numpy.amax(d["10"])) 
					
					block_ID = "block_" + str(start) #name the new block
					d = {} #empty dictionary
					GT = "" #start new GT
					af = [] #empty allele frequency records
					switches = 0
					 
				#this is the first block
				else:
					block_ID = "block_" + str(start) #name the first block
					GT = "" #start new GT
					switches = 0
			else:
				if "1/2" not in line:
					l = line.split("\t")
					gt = l[1] + l[2]
					if gt not in d.keys():
						d[gt] = []
						
					#this is the first line in the block		
					if GT == "":
						count = 1
						chromosome = l[3]
						start_coord = l[4]
						end_coord = str(int(l[4]) + int(length))
						
					#this is not the first line in the block
					else:
						if gt == GT: #my gt is similar to the previous one
							count += 1 #I need to increase my count
						else: #I am changing genotype
							d[GT].append(count) #I need to save in the previous genotype the total count
							count = 1 #and for my new genotype count returns to 1
							switches += 1
					GT = gt #GT is replaced by the last gt
					
					#get allele frequency
					field = l[7].split(":")[1]
					val = field.split(",")
					
					if float(val[0]) != 0 and float(val[1]) != 0:
						allele_freq = float(val[0])/(float(val[0])+float(val[1]))
					else:
						print "WARNING!!! Heterozygous SNP with an allele supported by 0 reads!"
						print line
						allele_freq = 0.5
						
					af.append(allele_freq)
				
	#checking the last block the last block
	d[GT].append(count) #I need to append the last count
	if "01" in d.keys() and "10" in d.keys():
		if len(d["01"]) > 0 and len(d["10"]) > 0: #both haplotypes were observed in the previous block
			if numpy.amax(d["01"]) >= int(min_var) and numpy.amax(d["10"]) >= int(min_var): #both genotypes present blocks with more than 5 variants phased
				if numpy.mean(af) < 0.6 and numpy.mean(af) > 0.4:
					print chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + block_ID + "\t" + tag + "\t" + str(numpy.amax(d["01"])) + "\t" + str(numpy.amax(d["10"])) #report the results of the previous block
					print >>out, chromosome + "\t" + start_coord + "\t" + end_coord + "\t" + block_ID + "\t" + tag + "\t" + str(numpy.amax(d["01"])) + "\t" + str(numpy.amax(d["10"])) 
				
	hapcut_file.close()
	out.close()
	print "\n"
							
parser = argparse.ArgumentParser(description="Detect crossing-over events using HapCut2")
parser.add_argument("-bam", "--bam", dest="bam", action="store", help="BAM file with read alignment")
parser.add_argument("-v", "--vcf", dest="vcf", action="store", help="VCF file with the HETEROZYGOUS variants to be considered")
parser.add_argument("-t", "--tag", dest="tag", action="store", help="Tag to use in output files")
parser.add_argument("-m", "--min_var", dest="min_var", action="store", default=10, help="Minumum number of followed phased variants to consider a block [10]")
parser.add_argument("-s", "--max_switch", dest="max_switch", action="store", default=0.01, help="Maximum number of genotype switches per phased positions to consider a block [0.01")
args = parser.parse_args()

cross_over(args.bam,args.vcf,args.tag,args.min_var,args.max_switch)
