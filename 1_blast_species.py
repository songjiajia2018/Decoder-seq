#!/usr/bin/env python
#coding=utf-8

import os
import re
import argparse

parser = argparse.ArgumentParser(description='''Blast checks fastq sequencing data for contamination and its source and percentage''')

parser.add_argument('-f', '--fa', type=str, required=True, help='fa file')

args = parser.parse_args()


fa = args.fa
prefix = re.sub('.fa','',fa)
#os.system('''zcat %s | head -n 20000 | awk '{if(NR%%4==1){print ">"$1}else if(NR%%4==2){print $0}}'|sed 's/@//g' > %s''' %(fastq_gz, fa))
os.system("blastn -query %s -out %s_blast_result.xml -max_target_seqs 1 -outfmt 5 -db /home/songjia/reference/NT/nt -num_threads 20 -evalue 1e-5" %(fa,prefix))



##Extract gi number from XML result file, including: 
#Iteration_query-def: reads id; 
#Hit_id: matching sequence's gi number; 
#Hit_def: matching sequence species.
from collections import defaultdict

xmlfile=open(prefix+"_blast_result.xml","r")
outfile=open(prefix+"_xml_extract_gi.txt","w")

dict1=defaultdict(list)
for lines in xmlfile:
	line=lines.strip()
	read_id = re.match('<Iteration_query-def>.*</Iteration_query-def>',line)
	Hit_id = re.match('<Hit_id>.*</Hit_id>',line)
	Hit_def = re.match('<Hit_def>.*</Hit_def>',line)
	
	if read_id !=None:
		read_id=read_id.group()
		read_id = read_id.split("<")[1].split(">")[1]
		key=read_id

	elif Hit_id !=None:
		Hit_id = Hit_id.group()
		Hit_id = Hit_id.split("<")[1].split(">")[1]
		dict1[key].append(Hit_id)

	elif Hit_def !=None:
		Hit_def = Hit_def.group()
		Hit_def = Hit_def.split("<")[1].split(">")[1]
		dict1[key].append(Hit_def)

for key in dict1:
	outfile.write(key + "\t" + "\t".join(dict1[key])+"\n")



##Scientific names. The blast result yields the gi number, taxid, and species scientific name.
extract_gi = open(prefix+"_xml_extract_gi.txt","r")
gi2taxid = open("/home/songjia/reference/Taxonomy/cutted_nucl_gb.accession2taxid","r")
taxid2name = open("/home/songjia/reference/Taxonomy/names.dmp","r")
get_name = open(prefix+"_species_name.txt","w")

taxid_name_dict={}
for lines in taxid2name:
	if "scientific name" in lines:
		line = lines.strip().split("|")
		taxid = line[0].strip()
		name = line[1].strip()
		taxid_name_dict[taxid]=name
		
extract_dict=defaultdict(list)
for lines in extract_gi:
	line = lines.strip().split("\t")
	gi = line[1].split("|")[1]
	extract_dict[gi].append("\t".join(line))


gi_taxid_dict={}
for lines in gi2taxid:
	line = lines.strip().split("\t")
	GI = line[1]
	taxid = line[0]
	gi_taxid_dict[GI]=taxid

jiaoji=set(extract_dict.keys())&set(gi_taxid_dict.keys())

tax_list=taxid_name_dict.keys()

extract_gi = open(prefix+"_xml_extract_gi.txt","r")
for lines in extract_gi:
	line = lines.strip().split("\t")
	gi = line[1].split("|")[1]
	if gi in jiaoji:
		taxid=gi_taxid_dict[gi]
		if taxid in tax_list:
			get_name.write("\t".join(line)+"\t"+taxid_name_dict[taxid]+"\n")



##Calculate the ratio of each species' reads to the total number of extracted reads
from collections import Counter

species_name=open(prefix+"_species_name.txt","r")
blast_final_result =open(prefix+"_blast_final_result.txt","w")

name_list_all=[]
for lines in species_name:
	line = lines.strip().split("\t")
	name = line[-1]
	name_list_all.append(name)

count_result = Counter(name_list_all)
count_list = count_result.items()
count_list = sorted(count_list,key=(lambda x:x[1]),reverse=True)

blast_final_result.write("Name\tHit_reads\tpercentage1\tpercentage2\n")
for i in count_list:
	name = i[0]
	number = i[1]
	reads_num = 5000
	percentage1 = "%.2f%%"%(100*float(number)/float(reads_num))
	percentage2 ="%.2f%%"%(100*float(number)/float(len(name_list_all)))
	blast_final_result.write(name+"\t"+str(number)+"\t"+str(percentage1)+"\t"+str(percentage2)+"\n")
