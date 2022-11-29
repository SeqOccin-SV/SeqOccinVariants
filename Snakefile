#!/usr/bin/env python


configfile: "config.yaml"


include : "rules/methods.smk"


# Get's the absolute path of config['samples'] to get around the workdir
config['samples'] = abs_path(config['samples'])
# path = workflow.basedir


samples = pd.read_table(config['samples'], comment='#').set_index('sample', drop=False)


workdir: config["workdir"]


rule all:
	input:
		in1 = expand("calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz.tbi", sample=samples.index, tech=config['datatype'], mapping = get_mapping(), tools=get_tools()),
		html = "stats/pipeline_output.html",
		htmls = expand("stats/{sample}-{tech}-{mapping}-{svtools}-{snptools}-stats.html", sample=samples.index, tech=config['datatype'], mapping = get_mapping(), svtools=get_sv_tool(), snptools=get_snp_tool())


include : "rules/mapping.smk"
include : "rules/sv_calling.smk"
include : "rules/snp_calling.smk"
include : "rules/vcf_handling.smk"
include : "rules/html_prep.smk"
include : "rules/html_builder.smk"
include : "rules/single_html.smk"
include : "rules/bed_building.smk"








