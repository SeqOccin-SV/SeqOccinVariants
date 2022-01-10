# Pacbio variants detection pipeline

Snakemake pipeline for full variant detection on PacBio Data (CLR, CCS)

## Prerequisite to run PacBio variants detection pipeline on your dataset

To run this snakemake pipeline, you need to have access to a conda environment with snakemake >= 5 installed.
On Genotoul, you can simply load the following module

```bash
module load bioinfo/snakemake-5.25.0
```

Do not forget to deactivate your own environment if you decide to use this solution

Detecting SNPs and Indels on HiFi datasets use DeepVariant and also require access to Singularity which is handle by the pipeline on genotoul cluster

## Running the pipeline

### Snakemake

First, modify the ```config.yaml```. There are several element that can be modified.

- **ref** path to the reference fasta file. Will need corresponding fasta index
- **samples** path to the tsv file following one of the two possible content (see below)
- **datatype** Either CLR or CCS. If you want to use the pipeline on ONT data, use CLR
- **max-length** max-length for insertion and duplication detection. Large value increase RAM usage. Between 10000 and 50000 is a good base.
- **annotation** path to a fasta file for SV annotation. Annotated variants will have a SVANN field with the matching fasta header.


Then, you can modify either the example ```samples.tsv``` or ```bamSamples.tsv``` depending on your needs.

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


Finally, use the launch script

```bash
sbatch -c 1 -p workq -J pbv --mail-type=END,FAIL ./launch.pbvariants.smkj
```

If some jobs are failing, you will probably have to modify the content of the ```cluster.yaml``` to ask for more ressources.
