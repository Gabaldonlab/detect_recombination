#!/usr/bin/env python

"""
By Veronica Mixao
@gabaldonlab, BSC, Barcelona
vmixao@gmail.com

Objective: Detect recombination events when performing read-mapping to a phased reference genome
Last change: 16/01/2020
"""

import argparse
import os

##############################
#Edit the path of each program
bedtools = "/usr/bin/bedtools"
samtools = "/usr/bin/samtools"
##############################



def bam2bed(bamfile,tag,multimap):
	#Obtain a bed file with read pairs mapped in different chromosomes
	print "\t", "Obtaining pairs of reads mapped in different chromosomes..."
	if multimap == "no":
		os.system("samtools view  -h -q 60 " + bamfile + " | grep -v -e 'XA:Z:' -e 'SA:Z:' | awk '($3!=$7 && $7!=\"=\")'> " + tag + ".sam") #exclude MQ < 60, multi-mappings, and reads with both pairs aligning in the same chromosome
	
	else:
		os.system("samtools view  -h " + bamfile + " | awk '($3!=$7 && $7!=\"=\")'> " + tag + ".sam") #exclude reads with both pairs aligning in the same chromosome
	
	
	print "\t", "Converting to BED format..."
	os.system("samtools view -bS " + tag + ".sam | bedtools bamtobed -i stdin -tag NM > " + tag + ".reads.bed") #change to bed format and include the number of mismatches
	
	os.system("rm " + tag + ".sam")
	print "\t", "Warning!!! SAM file removed to save space!"
	


def detect_recomb(readL,maxM,tag,codeA,codeB,snps,type_rec):
	special_reads = {}
	if snps != "":
		#identify reads overlapping high-density SNP regions
		os.system(bedtools + " intersect -a " + tag + ".reads.bed" + " -b " + snps + " -c  > " + tag + ".reads.snp_count.bed")
		f_open = open(tag + ".reads.snp_count.bed")
		f = f_open.readlines()
		
		count_snps = 0
			
		for line in f:
			lin = line.split("\n")
			l = line.split("\t")
			
			if int(l[6]) > 0:
				count_snps += 1
				special_reads[l[3]] = l[6]		
		f_open.close()
		print "\t" + "Number of reads overlapping SNPs: " + str(count_snps)
		
	#Get information for each read
	bed_file = open(tag + ".reads.bed", "r")
	bed = bed_file.readlines()
	
	read_coord = {}
	counter = 0
	
	#get read information and save it in dictionary
	print "\t", "Getting read information..."
	
	for line in bed:
		lin = line.split("\n")
		l = lin[0].split("\t")
		
		chrm = l[0]
		start = l[1]
		end = l[2]
		name = l[3].split("/")[0]
		read_length = float(float(end)-float(start))
		mismatches = l[4]
		strand = l[5]
		
		if read_length >= int(readL): #make sure that a certain number of bp aligned
			if int(mismatches) <= int(maxM): #make sure that they do not have many mismatches
				counter += 1
				if name not in read_coord.keys():
					read_coord[name] = []
				if codeA in chrm:
					species = "A"
				else:
					if codeB in chrm:
						species = "B"
					else:
						"Please revise the code of your chromosomes name!!!"
				info = chrm + "\t" + start + "\t" + end + "\t" + str(mismatches) + "\t" + name + "\t" + strand + "\t" + species
				read_coord[name].append(info)
			
			else:
				if l[3] in special_reads.keys():
					n_snps = int(special_reads[l[3]])
					if int(mismatches) == int(n_snps):
						counter += 1
						if name not in read_coord.keys():
							read_coord[name] = []
						if codeA in chrm:
							species = "A"
						else:
							if codeB in chrm:
								species = "B"
							else:
								"Please revise the code of your chromosomes name!!!"
						info = chrm + "\t" + start + "\t" + end + "\t" + str(mismatches) + "\t" + name + "\t" + strand + "\t" + species
						read_coord[name].append(info)
			
	bed_file.close()			
					
	print "\t", "\t", "Reads with potential interest: " + str(counter)
	
	if type_rec == "inter":
		#get read-pairs mapping in chromosomes of different species
		out = open(tag + ".recombination.reads.bed", "w+")
	
		i = 0
		for read in read_coord.keys():
			if len(read_coord[read]) == 2: #make sure that the two reads of a pair passed the filters
				chrm1,start1,end1,mismatches1,name1,strand1,species1 = read_coord[read][0].split("\t")
				chrm2,start2,end2,mismatches2,name2,strand2,species2 = read_coord[read][1].split("\t")
				
				if species1 != species2: #make sure that they align in chromosomes from different species				
					i += 1
					print >>out, chrm1 + "\t" + start1 + "\t" + end1 + "\t" + mismatches1 + "\t" + name1 + "\t" + strand1 + "\t" + chrm2 + "\t" + start2 + "\t" + end2 + "\t" + mismatches2 + "\t" + name2 + "\t" + strand2
					
		print "\t", "\t", "Reads supporting inter-species recombination: ", str(i)
		out.close()
	else:
		if type_rec == "intra":
			#get read-pairs mapping in chromosomes of the same species
			out = open(tag + ".rearrangements.reads.bed", "w+")
	
			i = 0	
			for read in read_coord.keys():
				if len(read_coord[read]) == 2: #make sure that the two reads of a pair passed the filters
					chrm1,start1,end1,mismatches1,name1,strand1,species1 = read_coord[read][0].split("\t")
					chrm2,start2,end2,mismatches2,name2,strand2,species2 = read_coord[read][1].split("\t")
				
					if species1 != species2: #make sure that they align in chromosomes from different species				
						i += 1
						print >>out, chrm1 + "\t" + start1 + "\t" + end1 + "\t" + mismatches1 + "\t" + name1 + "\t" + strand1 + "\t" + chrm2 + "\t" + start2 + "\t" + end2 + "\t" + mismatches2 + "\t" + name2 + "\t" + strand2
					
			print "\t", "\t", "Reads supporting inter-species recombination: ", str(i)
			out.close()
		else:
			print "Type of recombination is not valid!"


	
def filter_results(tag,d,minCov,type_rec):
	
	if type_rec == "inter":
		in_name = tag + ".recombination.reads.bed"
	if type_rec == "intra":
		in_name = tag + ".rearrangements.reads.bed"
	
	print "\t", "Getting recombination blocks and filtering the results..."
	os.system("sort -k1,1 -k 2,2n " + in_name + " > " + in_name + "Sorted") #sort bed file
	os.system(bedtools + " merge -c 5,6 -o count,collapse -d " + str(d) + " -i " + in_name + "Sorted > temp1.bed") #merge results for species 1
	os.system("rm " + in_name) #we do not need this file anymore
	
	os.system(bedtools + " intersect -a temp1.bed -b " + in_name + "Sorted -wo > temp2.bed") #recover information on the matches of each species1 block on species2
	os.system("rm temp1.bed")
	
	temp_open = open("temp2.bed", "r")
	temp = temp_open.readlines()
	temp3 = open("temp3.bed", "w+") 
	
	for line in temp:
		lin = line.split("\n")
		l = lin[0].split("\t")
		print >>temp3, l[11] + "\t" + l[12] + "\t" + l[13] + "\t" + l[16] + "\t" + l[0] + "\t" + l[1] + "\t" + l[2] + "\t" + l[3] + "\t" + l[4] #invert the order of the species so now we can merge the blocks of species2
	temp_open.close()
	os.system("rm temp2.bed")
	
	os.system("sort -k1,1 -k 2,2n temp3.bed > temp3.bedSorted")
	
	os.system(bedtools + " merge -c 4,4 -o count,collapse -d " + str(d) + " -i temp3.bedSorted > temp4.bed") #repeating but for species 2
	os.system(bedtools + " intersect -a temp4.bed -b temp3.bedSorted -wao > temp5.bed") 
	
	os.system("rm temp3.bed temp3.bedSorted temp4.bed")
	
	#Getting the final results
	simplified_results_file = open("temp5.bed", "r")
	simpl = simplified_results_file.readlines()
	final_out = open(tag + "." + type_rec + "_events.txt", "w+")
	final_out_bed = open(tag + "." + type_rec + "_events.bed", "w+")
	
	printed = []
	recomb = 0 #total number of events
	
	for line in simpl:
		lin = line.split("\n")
		l = lin[0].split("\t")
		chr1 = l[0]
		start1 = l[1]
		end1 = l[2]
		reads1 = l[3]
		strand1 = l[4]
		chr2 = l[9]
		start2 = l[10]
		end2 = l[11]
		reads2 = l[12]
		strand2 = l[13]
		
		counter_1plus = 0 #determine the strand of each alinment in species 1
		counter_2plus = 0 #determine the strand of each alinment in species 2
		
		if float(reads1) >= float(minCov) and float(reads2) >= float(minCov): #make sure that a given region has more than X reads aligning to a different species
			for strandness in strand1.split(","): 
				if strandness == "+":
					counter_1plus += 1
			if counter_1plus/len(strand1.split(",")) > 0.9:
				pos1 = end1 #if the majority of reads align in positive strand this means that the break is on the right
			else:
				if counter_1plus/len(strand1.split(",")) < 0.1:
					pos1 = start1 #if the majority of reads supporting the recombination are in the negative strand this means that the break is on the left
				else:
					pos1 = start1 + "-" + end1 #if we cannot stablish correctly the strand, the scripts outputs the block
			
			for strandness in strand2.split(","):
				if strandness == "+":
					counter_2plus += 1
			if counter_2plus/len(strand2.split(",")) > 0.9:
				pos2 = end2
			else:
				if counter_2plus/len(strand2.split(",")) < 0.1:
					pos2 = start2
				else:
					pos2 = start2 + "-" + end2
			
			info = chr1 + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2 + "\t" + chr2 + "\t" + start2 + "\t" + end2 + "\t" + chr1 + "\t" + start1 + "\t" + end1
			
			infoA = chr1 + "\t" + pos1 + "\t" + chr2 + "\t" + pos2 #final information
			info_bedA = chr1 + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2
			
			infoB = chr2 + "\t" + pos2 + "\t" + chr1 + "\t" + pos1 #final information
			info_bedB = chr2 + "\t" + start2 + "\t" + end2 + "\t" + chr1 + "\t" + start1 + "\t" + end1
			
			if info not in printed:
				print >>final_out, infoA
				print >>final_out, infoB
				print >>final_out_bed, info_bedA
				print >>final_out_bed, info_bedB
				printed.append(info)
				recomb += 1
		
	simplified_results_file.close()
	os.system("rm temp5.bed")
	print "\t" + "\t" + "\t" + "\t" + "TOTAL NUMBER OF RECOMBINATION EVENTS DETECTED FOR " + tag + ": " + str(recomb)
	print "DONE!!"
	
parser = argparse.ArgumentParser(description="Detect crossing-over events in samples mapped to a phased genome")
parser.add_argument("-bam", "--bam", dest="bam", action="store", help="BAM file with read alignment")
parser.add_argument("-l", "--length", dest="readL", action="store", default=75, help="Minimum number of bases for a read to be considered [75]")
parser.add_argument("-max", "--maxM", dest="maxM", action="store", default=2, help="Maximum number of mismatches allowed per read [2]")
parser.add_argument("-t", "--tag", dest="tag", action="store", default="strain",help="TAG for intermediate files [strain]")
parser.add_argument("-codeA", "--code_chrm_A", dest="codeA", action="store", help="Expression that can be used to identify chromosomes of species A, i.e. something that is not present in species B chromosome name, but is shared by all species A chromosome name")
parser.add_argument("-codeB", "--code_chrm_B", dest="codeB", action="store", help="Expression that can be used to identify chromosomes of species B, i.e. something that is not present in species A chromosome name, but is shared by all species B chromosome name")
parser.add_argument("-d", "--dist", dest="dist", action="store", default=100,help="Maximum distance allowed for two intervals be merged in the Summary file [100]")
parser.add_argument("-cov", "--coverage", dest="minCov", action="store", default=30, help="Minimum number of read pairs supporting the cross-over event [30]")
parser.add_argument("-mp", "--multi_map", dest="multimap", action="store", default="no", help="Consider reads with multiple alignments: no/yes [no]")
parser.add_argument("-v", "--vcf", dest="snps", action="store", default="", help="VCF file with SNPs that should be taken into account for calculation of mismatches")
parser.add_argument("-type", "--type_rec", dest="type_rec", action="store", default="inter", help="Type of recombination: 'inter' for inter-species or 'intra' for intra-species [inter]")
args = parser.parse_args()

print "\n******************************************************"
print "Preparing to detect recombination in ", args.tag, "\n"

bam2bed(args.bam,args.tag,args.multimap)
detect_recomb(args.readL,args.maxM,args.tag,args.codeA,args.codeB,args.snps,args.type_rec)
filter_results(args.tag,args.dist,args.minCov,args.type_rec)
