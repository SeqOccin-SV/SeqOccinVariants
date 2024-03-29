#!/usr/bin/env python


# ref = '/home/adf/seqoccin/data/reference/hs37d5_hsa10.fa'
import pandas as pd
from pprint import pprint
from pathlib import Path
import re
# import pathlib

###############################################################
### Methods
###############################################################


# def get_tools() :
# 	"""Get tools to use depending on the input data type.

# 	Return:
# 		list : List of tools to use depending on the data type.

# 	"""
# 	if config['datatype'] == 'CLR' :
# 		return ['longshot','pbsv']
# 	elif config['datatype'] == 'CCS' :
# 		return ['dv','pbsv']
# 	elif config['datatype'] == 'ONT' :
# 		return ['longshot','svim']
# 	return ['pbsv']



def get_sv_tool() :
	"""Get the tool used for sv calling depending on data type."""
	if config['datatype'] in ['CLR', 'CCS', 'hifi', 'ONT', 'nanopore'] :
		return 'cuteSV'
	return None


def get_snp_tool() :
	"""Get the tool used for snp calling depending on data type."""
	if config['datatype'] == 'CLR' :
		return 'longshot'
	elif config['datatype'] in ['CCS', 'hifi', 'ONT', 'nanopore'] :
		return 'dv'
	return None


def get_tools() :
	"""Get tools to use depending on the input data type.

	Return:
		list : List of tools to use depending on the data type.

	"""
	# return [get_snp_tool(), get_sv_tool()]
	return [get_snp_tool(), 'cuteSV']


def get_mapping() :
	"""Get tools to use depending on the input data type.

	Return:
		list: List of tools to use depending on the data type.

	"""
	if config['datatype'] in ['CLR', 'CCS'] :
		return 'pbmm2'
	elif config['datatype'] in ['ONT', 'nanopore'] :
		return 'minimap'
	return 'unknown'


def get_files(wildcards):
	"""Parse the samples.txt from config file.

	Args:
		wildcards (obj): object containing all of the snakemake wildcards.
	
	Return:
		list: list of files path.
	
	"""
	print(wildcards.sample)
	files = samples.loc[wildcards.sample, 'path'].split(',')
	return(files)


def get_bam(wildcards):
	"""Shortcut to use bam directly."""
	if 'bam_path' in samples.columns:
		return samples.loc[wildcards.sample, "bam_path"]
	#return "mapping/%s-%s-pbmm2.bam" % (wildcards.sample, wildcards.tech)
	return "mapping/%s-%s-%s.bam" % (wildcards.sample, wildcards.tech, wildcards.mapping)


def get_bai(wildcards):
	"""Shortcut to use bam directly."""
	if 'bam_path' in samples.columns:
		return samples.loc[wildcards.sample, "bam_path"]+'.bai'
	# return "mapping/%s-%s-pbmm2.bam.bai" % (wildcards.sample, wildcards.tech)
	return "mapping/%s-%s-%s.bam.bai" % (wildcards.sample, wildcards.tech, wildcards.mapping)


def get_suffix(string):
	"""Get the path file suffixes."""
	return(''.join(Path(string).suffixes))
	# return(''.join(pathlib.Path(string).suffixes))


def get_threads(rule, default):
	"""Get the default number of threads for a rule."""
	cluster_config = snakemake.workflow.cluster_config
	if rule in cluster_config and "threads" in cluster_config[rule]:
		return cluster_config[rule]["threads"]
	elif "default" in cluster_config and "threads" in cluster_config["default"]:
		return cluster_config["default"]["threads"]
	return default


def get_mem(rule, default):
	"""Get memory to allocate in Gigabytes."""
	cluster_config = snakemake.workflow.cluster_config
	if rule in cluster_config and "mem" in cluster_config[rule]:
		return cluster_config[rule]["mem"]
	elif "default" in cluster_config and "mem" in cluster_config["default"]:
		return cluster_config["default"]["mem"]
	return default


def get_mem_per_CPU(rule, default, threads=1):
	"""Calculate memory to allocate to each CPU

	The memory is converted from Gigabbytes to Megabytes to avoid type errors 
	(float/int) when called by samtools sort.
	"""
	cluster_config = snakemake.workflow.cluster_config
	if rule in cluster_config and "mem" in cluster_config[rule] and "threads" in cluster_config[rule] :
		return int((cluster_config[rule]["mem"]/cluster_config[rule]["threads"])*1000)
	elif "default" in cluster_config and "mem" in cluster_config["default"] and "threads" in cluster_config["default"] :
		return int((cluster_config["default"]["mem"]/cluster_config["default"]["threads"])*1000)
	return int((default/threads)*1000)


def get_fofn_types(fofn_file):
	"""Get the type of data (bam/fastq) of the files whose names are in the fofn file."""
	pprint('Inside get_fofn_types')
	type = 0
	with open(str(fofn_file), 'r') as fh:
		first_file = fh.readline().rstrip('\n')
		pprint('first_file : '+first_file)
		suffix = get_suffix(first_file)
		if (suffix.endswith('subreads.bam')):
			type = 1
		elif (suffix.endswith('ccs.bam')):
			type = 2
		elif (suffix.endswith('fastq.gz') or suffix.endswith('fq.gz') or suffix.endswith('fastq') or suffix.endswith('fq')):
			type = 3
	pprint('Type value: '+str(type))
	return(type)


def abs_path(file_name):
	"""Get absolute path to a file from it's relative path.
	
	this function has no effect if an absolute path is given.
	
	"""
	return str(Path(file_name).resolve())


def get_fai() :
	"""Return the path to the fai file associated with fasta ref."""
	return config['ref'] + '.fai'


def get_bed_name() :
	"""Use the referance file to name the associated bed file."""
	file = config['ref'].split('/')[-1]
	return '.'.join(file.split('.')[:-1]) + '.bed'


def get_fasta_name() :
	"""Get the name of the fasta file without the extention."""
	return '.'.join(config['ref'].split('/')[-1].split('.')[:-1])


def get_chromosomes() :
	is_chr = re.compile(r'(chr)?[0-9,X]+$')
	chrs = ['unmapped']
	with open(get_fai(), 'r') as fai :
		for line in fai :
			name = line.split()[0]
			if re.match(is_chr, name) :
				chrs.append(name)
	return chrs


def bed_split_name(bed, chr) :
	basename = '.'.join(bed.split('.')[:-1])
	return basename + "-" + chr + ".bed"



def split_bed() :
	bed = get_bed_name()
	chrs = get_chromosomes()
	with open(bed, 'r') as file :
		with open(bed_split_name(bed, "unmapped"), 'a') as unmapped :
			for line in file :
				chr = line.split()[0]
				if chr in chrs :
					with open(bed_split_name(bed, chr), 'w') as chr_bed :
						chr_bed.write(line)
				else :
					unmapped.write(line)



def get_cutesv_param() :
	"""Get param for cuteSV depending on sequencing tech following recomandations."""
	if config['datatype'] == 'CLR' :
		return '--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5'
	elif config['datatype'] == 'CCS' :
		return '--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5'
	elif config['datatype'] == 'ONT' :
		return '--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3'
	return 'unknown'















