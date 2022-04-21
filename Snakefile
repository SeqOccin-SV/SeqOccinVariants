#!/usr/bin/env python


configfile: "config.yaml"


include : "rules/methods.smk"


# Get's the absolute path of config['samples'] to get around the workdir
config['samples'] = abs_path(config['samples'])


samples = pd.read_table(config['samples'], comment='#').set_index('sample', drop=False)


workdir: config["workdir"]


rule all:
	input:
		expand("calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz.tbi", sample=samples.index, tech=config['datatype'], mapping = get_mapping(), tools=get_tools())


#rule test:
	#input:
		#bam = get_bam
	#output:
		#"calling/goat-ONT-minimap-svim.vcf.gz.tbi"
	#run:
		#print({input.bam})
		#f = open({output}, "a")
		#f.write("hello World")
		#f.close()

include : "rules/mapping.smk"
include : "rules/sv_calling.smk"
include : "rules/snp_calling.smk"
include : "rules/vcf_handling.smk"











