#!/usr/bin/env python


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


# Build stats files from cuteSV output files
rule cuteSV_build_stats:
    input:
        "calling/{sample}-{tech}-{mapping}-cuteSV.vcf.gz"
    output:
        "stats/{sample}-{tech}-{mapping}-cuteSV.vcf.stats.txt"
    conda:
        '../envs/html_env.yaml'
    shell:
        """
        bcftools query -f "%ID\n" {input} | cut -f2 -d'.' | sort | uniq -c > {output}
        """


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


# Index fofn file with corresponding sam.stats pbmm2 output file
rule index_fofn_bamstats:
    input:
        fofn = "fofn/{sample}-{tech}.fofn",
        bam = "mapping/{sample}-{tech}-{mapping}.bam.stats"
    output:
        temp("stats/{sample}-{tech}-{mapping}_index.txt")
    shell:
        "echo -e \"`pwd`/{input.bam}\\t{input.fofn}\" > {output}"



