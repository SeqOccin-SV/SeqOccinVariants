#!/usr/bin/env python


rule build_bed:
	input:
		fai = get_fai()
	output:
		bed = "{ref}_full.bed"
	shell:
		"""
		awk 'BEGIN {{FS = "\t"}}; {{print $1 FS \"0\" FS $2}}' {input.fai} > {output.bed}
		"""


rule trim_bed:
    wildcard_constraints:
        ref = get_fasta_name()
    input:
        bed = "{ref}_full.bed"
    output:
        trim_bed = "{ref}.bed"
    params:
        cutoff = config['cutoff']
    shell:
        """
        awk "{{ if (\$3 > {params.cutoff}) {{ print }} }}" {input.bed} > {output.trim_bed}
        """


rule split_bed:
    input:
        bed = expand("{ref}.bed", ref=get_fasta_name())
    output:
        chr_file = expand("{ref}-{chr}.bed", ref=get_fasta_name(), chr=get_chromosomes())
    run:
        split_bed()