#!/usr/bin/env python


rule single_bam_stats_for_html:
    input:
        "mapping/{sample}-{tech}-{mapping}.bam.stats"
    output:
        "stats/{sample}-{tech}-{mapping}.bam.html.stats"
    log:
        "logs/html/{sample}-{tech}-{mapping}.bam.html.stats.log"
    shell:
        "echo {input} > {output}"


# Groupe pbsv stats output
rule single_pbsv_stats_for_html:
    input:
        "stats/{sample}-{tech}-pbmm2-pbsv.vcf.stats.txt"
    output:
        # "stats/pbsv_stats_for_html.stats"
        "stats/{sample}-{tech}-pbmm2-pbsv.vcf.html.stats.txt"
    shell:
        "echo {input} > {output}"


# Groupe svim stats output
rule single_svim_stats_for_html:
    input:
        "stats/{sample}-{tech}-minimap-svim.vcf.stats.txt"
    output:
        # "stats/svim_stats_for_html.stats"
        "stats/{sample}-{tech}-minimap-svim.vcf.html.stats.txt"
    shell:
        "echo {input} > {output}"


# Groupe cuteSV stats output
rule single_cuteSV_stats_for_html:
    input:
        "stats/{sample}-{tech}-{mapping}-cuteSV.vcf.stats.txt"
    output:
        # "stats/svim_stats_for_html.stats"
        "stats/{sample}-{tech}-{mapping}-cuteSV.vcf.html.stats.txt"
    shell:
        "echo {input} > {output}"


# Groupe all variant sizes information
rule single_variantsizes_for_html:
    input:
        "stats/{sample}-{tech}-{mapping}-{svtools}.variantsizes.tsv"
    output:
        # "stats/{tools}_txt_variantsizes.tsv"
        "stats/{sample}-{tech}-{mapping}-{svtools}.html.variantsizes.tsv"
    shell:
        "echo {input} > {output}"


# Build the stats file for the html output from the outputs of longshot or deepvariant (CLR or CCS)        
if "longshot" in get_tools() :

    rule single_longshot_stats_for_html:
        input:
            "stats/{sample}-{tech}-{mapping}-longshot.vcf.gz.stats"
        output:
            "stats/{sample}-{tech}-{mapping}-longshot.vcf.gz.html.stats"
        shell:
            "echo {input} > {output}"
else :

    rule single_dv_stats_for_html:
        input:
            "stats/{sample}-{tech}-{mapping}-dv.vcf.gz.stats"
        output:
            # "stats/SNP_dv.stats"
            "stats/{sample}-{tech}-{mapping}-dv.vcf.gz.html.stats"
        shell:
            "echo {input} > {output}"


rule single_group_sequence_numbers:
    input:
        "stats/{sample}-{tech}_sequences_count.stats"
    output:
        # "stats/sequences_count.stats"
        "stats/{sample}-{tech}_sequences_count.html.stats"
    shell:
        "cat {input} > {output}"


rule single_group_index_fofn_bamstats:
    input:
        "stats/{sample}-{tech}-{mapping}_index.txt"
    output:
        "stats/{sample}-{tech}-{mapping}_index.html.txt"
    shell:
        "cat {input} > {output}"


# rule that calls the rules above to build the html output
rule build_single_html:
    input:
        bam = "stats/{sample}-{tech}-{mapping}.bam.html.stats",
        sv = "stats/{sample}-{tech}-{mapping}-{svtools}.vcf.html.stats.txt",
        variantsizes = "stats/{sample}-{tech}-{mapping}-{svtools}.html.variantsizes.tsv",
        snp = "stats/{sample}-{tech}-{mapping}-{snptools}.vcf.gz.html.stats",
        index = "stats/{sample}-{tech}-{mapping}_index.html.txt",
        number = "stats/{sample}-{tech}_sequences_count.html.stats"
    output:
        "stats/{sample}-{tech}-{mapping}-{svtools}-{snptools}-stats.html"
    log:
        "logs/stats/{sample}-{tech}-{mapping}-{svtools}-{snptools}-stats.html.log"
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

