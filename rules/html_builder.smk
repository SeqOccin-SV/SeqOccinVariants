#!/usr/bin/env python

# Groupe all outputs of samtools stats from the pbmm2 output files
rule bam_stats_for_html:
	input:
		expand("mapping/{sample}-{tech}-{mapping}.bam.stats", sample=samples.index, tech=config['datatype'], mapping = get_mapping())
	output:
		"stats/bam_stats_for_html.stats"
	log:
		"logs/html/bam_stats_for_html.stats.log"
	threads: 1
	shell:
		"echo {input} | tr \" \" \"\\n\" > {output}"



# Groupe pbsv stats output
rule pbsv_stats_for_html:
    input:
        expand("stats/{sample}-{tech}-pbmm2-pbsv.vcf.stats.txt", sample=samples.index, tech=config['datatype'])
    output:
        "stats/pbsv_stats_for_html.stats"
    shell:
        "echo {input} | tr \" \" \"\\n\" > {output}"


# Groupe svim stats output
rule svim_stats_for_html:
    input:
        expand("stats/{sample}-{tech}-minimap-svim.vcf.stats.txt", sample=samples.index, tech=config['datatype'])
    output:
        "stats/svim_stats_for_html.stats"
    shell:
        "echo {input} | tr \" \" \"\\n\" > {output}"


# Groupe all variant sizes information
rule variantsizes_for_html:
    input:
        expand("stats/{sample}-{tech}-{mapping}-{tools}.variantsizes.tsv", sample=samples.index, tech=config['datatype'], mapping = get_mapping(), tools=get_sv_tool())
    output:
        "stats/{tools}_txt_variantsizes.tsv"
    shell:
        "echo {input} | tr \" \" \"\\n\" > {output}"


# Build the stats file for the html output from the outputs of longshot or deepvariant (CLR or CCS)        
if "longshot" in get_tools() :

    rule longshot_stats_for_html:
        input:
            expand("stats/{sample}-{tech}-{mapping}-longshot.vcf.gz.stats", sample=samples.index, tech=config['datatype'], mapping = get_mapping())
        output:
            "stats/SNP_longshot.stats"
        shell:
            "echo {input} | tr \" \" \"\\n\" > {output}"
else :

    rule dv_stats_for_html:
        input:
            expand("stats/{sample}-{tech}-{mapping}-dv.vcf.gz.stats", sample=samples.index, tech=config['datatype'], mapping = get_mapping())
        output:
            "stats/SNP_dv.stats"
        shell:
            "echo {input} | tr \" \" \"\\n\" > {output}"



rule group_sequence_numbers:
    input:
        expand("stats/{sample}-{tech}_sequences_count.stats", sample=samples.index, tech=config['datatype'])
    output:
        "stats/sequences_count.stats"
    shell:
        "for i in {input}; do cat $i >> {output}; done"


rule group_index_fofn_bamstats:
    input:
        expand("stats/{sample}-{tech}-{mapping}_index.txt", sample=samples.index, tech=config['datatype'], mapping = get_mapping())
    output:
        "stats/fofn_bamstats_index.txt"
    shell:
        "cat {input} > {output}"


# rule that calls the rules above to build the html output
rule build_full_html:
    input:
        bam = "stats/bam_stats_for_html.stats",
        sv = expand("stats/{tools}_stats_for_html.stats", tools=get_sv_tool()),
        variantsizes = expand("stats/{tools}_txt_variantsizes.tsv", tools=get_sv_tool()),
        snp = expand("stats/SNP_{tools}.stats", tools=get_snp_tool()),
        index = "stats/fofn_bamstats_index.txt",
        number = "stats/sequences_count.stats"
    output:
        "stats/pipeline_output.html"
    log:
        "logs/stats/pipeline_output.html.log"
    conda:
        '../envs/html_env.yaml'
    params:
        path = workflow.basedir,
        dtype = config['datatype']
    shell:
        """
        python3 {params.path}/scripts/html_methods.py \
        -b {input.bam} \
        -p {input.sv} \
        -v {input.variantsizes} \
        -s {input.snp} \
        -i {input.index} \
        -n {input.number} \
        -t {output} \
        -w {params.path}/envs \
        -d {params.dtype} 2> {log}
        """
