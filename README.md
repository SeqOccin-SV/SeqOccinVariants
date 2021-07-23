# Pacbio variants detection pipeline

Snakemake pipeline for full variant detection on PacBio Data (CLR, CCS)

## Prerequisite to run PacBio variants detection pipeline on your dataset

To run this snakemake pipeline, you need to have access to a conda environment with snakemake >= 5 installed.
On Genotoul, you can simply load the following module

```bash
module load bioinfo/snakemake-5.25.0
```

Detecting SNPs and Indels on HiFi datasets use DeepVariant and also require access to Singularity

## Running the pipeline

### Snakemake

First, modify the ```config.yaml``` then either the example ```samples.tsv``` or ```bamSamples.tsv```

You can write your own tsv file from scratch if you respect the following content

Example samples.tsv
|sample|path|
|-------|----|
|name_CLR|1.subreads.bam,2.subreads.bam|
|name_HiFi1|1.ccs.bam,1.ccs.bam|
|name_HiFi2|1.fq.gz,2.fq.gz|

Example bamSamples.tsv

|sample|bam_path|
|-------|----|
|name|path/mapped.bam|

then use the launch script

```bash
sbatch -j 4 ./launch.pbsv.smkj
```
