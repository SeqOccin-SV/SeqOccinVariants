#!/usr/bin/env python


### VCF Handling

# calling/{sample}-{tech}-{mapping}-{tools}.vcf

ruleorder: compress_tabix > tabix

rule compress_tabix:
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
