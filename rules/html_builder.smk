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


# Build stats files from pbsv output files
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

# Groupe pbsv stats output
rule pbsv_stats_for_html:
    input:
        expand("stats/{sample}-{tech}-pbmm2-pbsv.vcf.stats.txt", sample=samples.index, tech=config['datatype'])
    output:
        "stats/pbsv_stats_for_html.stats"
    shell:
        "echo {input} | tr \" \" \"\\n\" > {output}"


# Build stats files from svim output files
rule svim_build_stats:
    input:
        "calling/{sample}-{tech}-minimap-svim.vcf.gz"
    output:
        "stats/{sample}-{tech}-minimap-svim.vcf.stats.txt"
    conda:
        '../envs/html_env.yaml'
    shell:
        """
        bcftools query -f "%ID\n" {input} | cut -f2 -d'.' | sort | uniq -c > {output}
        """

# Groupe svim stats output
rule svim_stats_for_html:
    input:
        expand("stats/{sample}-{tech}-minimap-svim.vcf.stats.txt", sample=samples.index, tech=config['datatype'])
    output:
        "stats/svim_stats_for_html.stats"
    shell:
        "echo {input} | tr \" \" \"\\n\" > {output}"


# Calculate variant sizes for each individual
rule set_variantsizes:
    input:
        "calling/{sample}-{tech}-{mapping}-{tools}.vcf.gz"
    output:
        "stats/{sample}-{tech}-{mapping}-{tools}.variantsizes.tsv"
    conda:
        '../envs/html_env.yaml'
    shell:
        """
        printf \"chrom\tpos\tsvtype\tsize\n\" > {output}; \
        bcftools query -i 'INFO/SVTYPE=="DEL" || INFO/SVTYPE=="INS"' -f "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n" {input} \
        >> {output}
        """


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



# # Get a file containing a list of impout files to pbmm2 with the number of sequences for each file.
# rule get_sequence_numbers:
#     input:
#         "fofn/{sample}-{tech}.fofn"
#     output:
#         temp("stats/{sample}-{tech}_sequences_count.stats")
#     conda:
#         '../envs/samtools_env.yaml'
#     shell:
#         """
#         for i in `cat {input}`
#         do 
#             sufix=`head -1 {input} | rev | cut -d'.' -f1 | rev`
#             if [ $sufix = \"bam\" ]
#             then
#                 echo -e \"{input}\\t`samtools view $i | wc -l`\" >> {output}
#             else
#                 echo -e \"{input}\\t`zcat $i | grep -c \"@\"`\" >> {output}
#             fi
#         done
#         """



# Get a file containing a list of impout files to pbmm2 with the number of sequences for each file.
rule get_sequence_numbers:
    input:
        "fofn/{sample}-{tech}.fofn"
    output:
        temp("stats/{sample}-{tech}_sequences_count.stats")
    conda:
        '../envs/samtools_env.yaml'
    shell:
        """
        for i in `cat {input}`
        do 
            sufix=`head -1 {input} | rev | cut -d'.' -f1 | rev`
            if [ $sufix = \"bam\" ]
            then
                echo -e \"{input}\\t`samtools stats $i | \
                sed -n '/\\tsequences:\\|\\ttotal length:/p' | \
                cut -f3 | tr \"\\n\" \" \" | sed -e 's/ \+/\\t/g'`\" >> {output}
            else
                echo -e \"{input}\\t`seqkit stats $i | \
                tail -1 | sed -e 's/ \\+ /\\t/g' | \
                cut -f4,5 | sed -e 's/,//g'`\" >> {output}
            fi
        done
        """


rule group_sequence_numbers:
    input:
        expand("stats/{sample}-{tech}_sequences_count.stats", sample=samples.index, tech=config['datatype'])
    output:
        "stats/sequences_count.stats"
    shell:
        "for i in {input}; do cat $i >> {output}; done"


# Index fofn file with corresponding sam.stats pbmm2 output file
rule index_fofn_bamstats:
    input:
        fofn = "fofn/{sample}-{tech}.fofn",
        bam = "mapping/{sample}-{tech}-{mapping}.bam.stats"
    output:
        temp("stats/{sample}-{tech}-{mapping}_index.txt")
    shell:
        "echo -e \"`pwd`/{input.bam}\\t{input.fofn}\" > {output}"


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
