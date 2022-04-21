#!/usr/bin/env python


### SNP Calling


# SNP calling for CLR data 
rule longshot:
	input:
		bam = get_bam,
		bai = get_bai,
		ref = config['ref'],
		stats = "mapping/{sample}-{tech}-pbmm2.bam.stats"
	wildcard_constraints:
		mapping = "pbmm2"
	output:
		"calling/{sample}-{tech}-{mapping}-longshot.vcf"
	log:
		"logs/longshot/{sample}-{tech}-{mapping}.log"
	conda:
		'../envs/longshot_env.yaml'
	threads: get_threads('pbsv_discover',10)
	shell:
		"""
		longshot --ref {input.ref} \
		--bam {input.bam} \
		--out {output} \
		2> {log}
		"""


#SNP calling for HiFi data
rule deepvariant:
	input:
		bam = get_bam,
		bai = get_bai,
		ref = config['ref'],
		stats = "mapping/{sample}-{tech}-pbmm2.bam.stats"
	wildcard_constraints:
		mapping = "pbmm2"
	output:
		vcf = "calling/{sample}-{tech}-{mapping}-dv.vcf.gz",
		gvcf = "calling/{sample}-{tech}-{mapping}-dv.gvcf.gz"
	log:
		"logs/deepvariant/{sample}-{tech}-{mapping}.log"
	# ~ container:
		# ~ 'docker://google/deepvariant:1.0.0'
	envmodules:
		"system/singularity-3.5.3"
	threads: get_threads('deepvariant',20)
	shell:
		"""
		singularity exec -B /usr/bin/locale:/usr/bin/locale,/work/project/seqoccin:/work/project/seqoccin \
		-W ./ docker://google/deepvariant:1.3.0 \
		/opt/deepvariant/bin/run_deepvariant --model_type=PACBIO \
		--ref={input.ref} \
		--reads={input.bam} \
		--output_vcf={output.vcf} \
		--output_gvcf={output.gvcf} \
		--sample_name={wildcards.sample} \
		--num_shards={threads}
		"""

