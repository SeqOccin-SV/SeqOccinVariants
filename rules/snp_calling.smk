#!/usr/bin/env python


### SNP Calling


# SNP calling for CLR data 
rule longshot:
	wildcard_constraints:
		tools = 'longshot'
	# 	mapping = 'pbmm2'
	input:
		bam = get_bam,
		bai = get_bai,
		ref = config['ref'],
		stats = "mapping/{sample}-{tech}-{mapping}.bam.stats"
	output:
		"calling/{sample}-{tech}-{mapping}-{tools}.vcf"
	log:
		"logs/longshot/{sample}-{tech}-{mapping}-{tools}.log"
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
	# wildcard_constraints:
	# 	mapping = "pbmm2"
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
		--num_shards={threads} \
        --intermediate_results_dir tmp/
		"""


# #PEPPER margin calling for ONT data
# rule pepper:
# 	input:
# 		bam = get_bam,
# 		bai = get_bai,
# 		ref = config['ref'],
# 		stats = "mapping/{sample}-{tech}-minimap.bam.stats"
# 	wildcard_constraints:
# 		mapping = "minimap"
# 	output:
# 		vcf = "calling/{sample}-{tech}-{mapping}-pepper.vcf.gz"
# 	params:
# 		prefix = "calling/{sample}-{tech}-{mapping}-pepper"
# 	log:
# 		"logs/pepper/{sample}-{tech}-{mapping}.log"
# 	envmodules:
# 		"system/singularity-3.5.3"
# 	threads: get_threads('deepvariant',20)
# 	shell:
# 		"""
# 		singularity exec --bind /usr/lib/locale:/usr/bin/locale,/work/project/seqoccin:/work/project/seqoccin \
# 		-W calling/ docker://kishwars/pepper_deepvariant:r0.8 \
# 		run_pepper_margin_deepvariant call_variant \
# 		-b {input.bam} \
# 		-f {input.ref} \
# 		-o calling \
# 		-p {params.prefix} \
# 		-t {threads} \
# 		--ont_r10_q20
# 		"""