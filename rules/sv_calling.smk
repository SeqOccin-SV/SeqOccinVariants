#!/usr/bin/env python


### SV CALLING


# first rule for sv detection, use bam to look for regions with possible variants
rule pbsv_discover:
	input:
		bam = get_bam,
		bai = get_bai,
		stats = "mapping/{sample}-{tech}-pbmm2.bam.stats"
	wildcard_constraints:
		mapping = "pbmm2"
	output:
		"calling/{sample}-{tech}-{mapping}.svsig.gz"
	log:
		"logs/pbsv/{sample}-{tech}_{mapping}-discover.log"
	conda:
		'../envs/pbsv_env.yaml'
	threads: get_threads('pbsv_discover',10)
	shell:
		"pbsv discover -s {wildcards.sample} {input.bam} {output} --log-file {log}"


# second rule for sv detection, use .svsig.gz information to call variants
rule pbsv_call:
	input:
		gz = "calling/{sample}-{tech}-pbmm2.svsig.gz",
		ref = config['ref']
	output:
		"calling/{sample}-{tech}-pbmm2-pbsv.vcf"
	resources:
		#type = rules.pbmm2.params.type
	params:
		#is_ccs = lambda wildcards, resources: '--ccs' if resources.type==2 or resources.type==3 else '',
		has_annotation = lambda wildcards, resources: '--annotations '+config['annotation'] if 'annotation' in config else '',
		max_length = str(config['max-length'])
	log:
		"logs/pbsv/{sample}-{tech}_call.log"
	benchmark:
		"bench/{sample}-{tech}.pbsvCall.benchmark.txt"
	conda:
		'../envs/pbsv_env.yaml'
	threads: get_threads('pbsv_call',20)
	shell:
		"""
		pbsv call -j {threads} \
		--max-ins-length {params.max_length} \
		--max-dup-length {params.max_length} \
		{params.has_annotation} {input.ref} {input.gz} {output} \
		2> {log}
		"""
		# specific parameter for CCS, easier to consider a SV" --ccs"
		# will not implemented it for now as it is usefull mainly for low cov


rule svim:
	input:
		#"mapping/{sample}-{tech}-minimap.bam"
		get_bam
	wildcard_constraints:
		mapping = "minimap"
	output:
		#"calling/{sample}-{tech}-minimap-{tools}.vcf"
		"calling/{sample}-{tech}-{mapping}-{tools}.vcf"
	params:
		genome = config['ref'],
		outdir = "calling"
	conda :
		'../envs/svim_env.yaml'
	shell:
		"svim alignment {params.outdir} {input} {params.genome};"
		"mv calling/variants.vcf {output};"
		"echo {input}"





















