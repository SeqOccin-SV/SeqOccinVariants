#!/usr/bin/env python


# ref = '/home/adf/seqoccin/data/reference/hs37d5_hsa10.fa'
import pandas as pd
from pprint import pprint
from pathlib import Path
# import pathlib

###############################################################
### Methods
###############################################################


#------------------------------------------------------------------------------
# info : Get tools to use depending on the input data type.
# return : [str] : List of tools to use depending on the data type.
#------------------------------------------------------------------------------
def get_tools() :
	if config['datatype'] == 'CLR' :
		return ['longshot','pbsv']
	elif config['datatype'] == 'CCS' :
		return ['dv','pbsv']
	elif config['datatype'] == 'ONT' :
		return ['svim']
	return ['pbsv']


#------------------------------------------------------------------------------
# info : Get tools to use depending on the input data type.
# return : [str] : List of tools to use depending on the data type.
#------------------------------------------------------------------------------
def get_mapping() :
	if config['datatype'] in ['CLR', 'CCS'] :
		return 'pbmm2'
	elif config['datatype'] in ['ONT', 'nanopore'] :
		return 'minimap'
	return 'unknown'


#------------------------------------------------------------------------------
# info : Parse the samples.txt from config file.
# input : obj : object containing all of the snakemake wildcards.
# return : [str] : list of files path.
#------------------------------------------------------------------------------
def get_files(wildcards):
	print(wildcards.sample)
	files = samples.loc[wildcards.sample, 'path'].split(',')
	return(files)


# Shortcut to use bam directly
def get_bam(wildcards):
	if 'bam_path' in samples.columns:
		return samples.loc[wildcards.sample, "bam_path"]
	#return "mapping/%s-%s-pbmm2.bam" % (wildcards.sample, wildcards.tech)
	return "mapping/%s-%s-%s.bam" % (wildcards.sample, wildcards.tech, wildcards.mapping)


# Shortcut to use bam directly
def get_bai(wildcards):
	if 'bam_path' in samples.columns:
		return samples.loc[wildcards.sample, "bam_path"]+'.bai'
	return "mapping/%s-%s-pbmm2.bam.bai" % (wildcards.sample, wildcards.tech)


# Return path file suffixes
def get_suffix(string):
	return(''.join(Path(string).suffixes))
	# return(''.join(pathlib.Path(string).suffixes))


def get_threads(rule, default):
	cluster_config = snakemake.workflow.cluster_config
	if rule in cluster_config and "threads" in cluster_config[rule]:
		return cluster_config[rule]["threads"]
	elif "default" in cluster_config and "threads" in cluster_config["default"]:
		return cluster_config["default"]["threads"]
	return default


#------------------------------------------------------------------------------
# info : Returning memory to allocate in Gigabytes.
# in : rule (str) : corresponding rule in cluster.yaml.
# in : default (int) : default amount of total memory in gigabytes.
# return : (int) : amount of memory to allocate in Gigabytes.
#------------------------------------------------------------------------------
def get_mem(rule, default):
	cluster_config = snakemake.workflow.cluster_config
	if rule in cluster_config and "mem" in cluster_config[rule]:
		return cluster_config[rule]["mem"]
	elif "default" in cluster_config and "mem" in cluster_config["default"]:
		return cluster_config["default"]["mem"]
	return default


#------------------------------------------------------------------------------
# info : Calculate memory to allocate to each CPU, converting from gigabytes 
# to megabytes to avoid type errors (float/int) when called by samtools sort.
# in : rule (str) : corresponding rule in cluster.yaml.
# in : default (int) : default amount of total memory in gigabytes.
# in : threads (int) : default number of threads.
# return : (int) : amount of memory per CPU in Megabytes.
#------------------------------------------------------------------------------
def get_mem_per_CPU(rule, default, threads=1):
	cluster_config = snakemake.workflow.cluster_config
	if rule in cluster_config and "mem" in cluster_config[rule] and "threads" in cluster_config[rule] :
		return int((cluster_config[rule]["mem"]/cluster_config[rule]["threads"])*1000)
	elif "default" in cluster_config and "mem" in cluster_config["default"] and "threads" in cluster_config["default"] :
		return int((cluster_config["default"]["mem"]/cluster_config["default"]["threads"])*1000)
	return int((default/threads)*1000)


def get_fofn_types(fofn_file):
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


#------------------------------------------------------------------------------
# Get absolute path to a file from it's relative path.
# This function has no effect if an absolute path is given.
#------------------------------------------------------------------------------
def abs_path(file_name):
	return str(Path(file_name).resolve())


















