#!/usr/bin/env python


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


rule pbsv_build_stats:
    input:
        "calling/{sample}-{tech}-pbmm2-pbsv.vcf.gz"
    output:
        "stats/{sample}-{tech}-pbmm2-pbsv.vcf.stats.txt"
    conda:
        '../envs/html_env.yaml'
    shell:
        """
        bcftools query -f "%ID\n" {input} | cut -f2 -d'.' | sort | uniq -c > {output}
        """

rule pbsv_stats_for_html:
    input:
        expand("stats/{sample}-{tech}-pbmm2-pbsv.vcf.stats.txt", sample=samples.index, tech=config['datatype'])
    output:
        "stats/pbsv_stats_for_html.stats"
    shell:
        "echo {input} | tr \" \" \"\\n\" > {output}"


rule set_variantsizes:
    input:
        "calling/{sample}-{tech}-{mapping}-pbsv.vcf.gz"
    output:
        "stats/{sample}-{tech}-{mapping}-pbsv.variantsizes.tsv"
    conda:
        '../envs/html_env.yaml'
    shell:
        """
        printf \"chrom\tpos\tsvtype\tsize\n\" > {output}; \
        bcftools query -i 'INFO/SVTYPE=="DEL" || INFO/SVTYPE=="INS"' -f "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n" {input} \
        >> {output}
        """


rule variantsizes_for_html:
    input:
        expand("stats/{sample}-{tech}-{mapping}-pbsv.variantsizes.tsv", sample=samples.index, tech=config['datatype'], mapping = get_mapping(), tools=get_tools())
    output:
        "stats/txt_variantsizes.tsv"
    shell:
        "echo {input} | tr \" \" \"\\n\" > {output}"

        
if "longshot" in get_tools() :
    rule bcf_stats_longshot:
        input:
            # expand("calling/{sample}-{tech}-{mapping}-longshot.vcf.gz", sample=samples.index, tech=config['datatype'], mapping = get_mapping())
            "calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz"
        output:
            "stats/{sample}-{tech}-{mapping}-{tools}.vcf.gz.stats"
        conda:
            '../envs/html_env.yaml'
        shell:
            "bcftools stats {input} > {output}"

    rule longshot_stats_for_html:
        input:
            expand("stats/{sample}-{tech}-{mapping}-longshot.vcf.gz.stats", sample=samples.index, tech=config['datatype'], mapping = get_mapping())
        output:
            "stats/SNP_longshot.stats"
        shell:
            "echo {input} | tr \" \" \"\\n\" > {output}"
else :
    rule bcf_stats_dv:
        input:
            # expand("calling/{sample}-{tech}-{mapping}-dv.vcf.gz", sample=samples.index, tech=config['datatype'], mapping = get_mapping(), tools=get_tools())
            "calling/{sample}-{tech}-{mapping}-dv.vcf.gz"
        output:
            "stats/{sample}-{tech}-{mapping}-{tools}.vcf.gz.stats"
        conda:
            '../envs/html_env.yaml'
        shell:
            "bcftools stats {input} > {output}"

    rule dv_stats_for_html:
        input:
            expand("stats/{sample}-{tech}-{mapping}-dv.vcf.gz.stats", sample=samples.index, tech=config['datatype'], mapping = get_mapping())
        output:
            "stats/SNP_dv.stats"
        shell:
            "echo {input} | tr \" \" \"\\n\" > {output}"


rule build_full_html:
    input:
        bam = "stats/bam_stats_for_html.stats",
        pbsv = "stats/pbsv_stats_for_html.stats",
        variantsizes = "stats/txt_variantsizes.tsv",
        snp = expand("stats/SNP_{tools}.stats", tools=get_tools()[0])
    output:
        "stats/pipeline_output.html"
    # run:
        # write_full_html(str(input.bam), str(input.pbsv), str(input.variantsizes), str(input.snp), html = str(output))
    log:
        "logs/stats/pipeline_output.html.log"
    conda:
        '../envs/html_env.yaml'
    params:
        path = workflow.basedir
    shell:
        "python3 {params.path}/scripts/html_methods.py -b {input.bam} -p {input.pbsv} -v {input.variantsizes} -s {input.snp} -t {output} 2> {log}"
