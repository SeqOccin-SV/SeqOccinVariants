#!/usr/bin/env python


### VCF Handling


##########################################################################
# bcftools sort with Longshot vcf is not consisten and can produce errors
# but simply applying a bgzip can leed to tabix making errors if the vcf
# file is not properly sorted, as observed with tools such as svim.
# To manage these problems bgzip is used when longshot is used and
# bcftools sort is used otherwise.
##########################################################################


ruleorder: compress_tabix_longshot > compress_tabix > tabix

rule compress_tabix_longshot:
	wildcard_constraints:
		tools = 'longshot'
	input:
		vcf = "calling/{sample}-{tech}-{mapping}-{tools}.vcf"
	output:
		zip = "calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz",
		idx = "calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz.tbi"
	log:
		"logs/tabix/{sample}-{tech}-{mapping}-{tools}.log"
	conda:
		'../envs/tabix_env.yaml'
	threads: 1
	shell:
		"bgzip {input.vcf}; "
		"tabix {output.zip}"


rule compress_tabix:
	wildcard_constraints:
		tools = '|'.join([x for x in get_tools() if x != 'longshot'])
	input:
		vcf = "calling/{sample}-{tech}-{mapping}-{tools}.vcf"
	output:
		zip = "calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz",
		idx = "calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz.tbi"
	log:
		"logs/tabix/{sample}-{tech}-{mapping}-{tools}.log"
	conda:
		'../envs/tabix_env.yaml'
	threads: 1
	shell:
		"bcftools sort -Oz {input.vcf} -o {output.zip}; "
		"tabix {output.zip}"


rule tabix:
	input:
		gz = "calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz"
	output:
		idx = "calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz.tbi"
	log:
		"logs/tabix/{sample}-{tech}-{mapping}-{tools}.log"
	conda:
		'../envs/tabix_env.yaml'
	threads: 1
	shell:
		"tabix {input}"
