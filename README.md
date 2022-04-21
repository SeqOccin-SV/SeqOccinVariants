# Pacbio variants detection pipeline

Snakemake pipeline for full variant detection on PacBio Data (CLR, CCS)

## Prerequisite to run PacBio variants detection pipeline on your dataset

To run this snakemake pipeline, you need to have access to a conda environment with snakemake >= 5 installed.
On Genotoul, you can simply load the following module :

```bash
module load bioinfo/snakemake-5.25.0
```

Do not forget to deactivate your own environment if you decide to use this solution.

Update : the module is loaded by default by the launcher, so simply executing the launcher with no other modules loaded on the cluster is the ideal solution.

Detecting SNPs and Indels on HiFi datasets use DeepVariant and also require access to Singularity which is handle by the pipeline on the genotoul cluster.

## Running the pipeline

### Snakemake

First, modify the ```config.yaml```. There are several element that can be modified.

- **ref** path to the reference fasta file. Will need corresponding fasta index (recommended for fasta indexation : samtools faidx)
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

Note : the tsv file can't contain mixed datatypes, having fastq and bam together, or CCS bam and CLR bam together, is not tolerated by pbmm2 and, by extension, by the pipeline.



Finally, use the launch script, on the genotoul cluster :

```bash
sbatch -c 1 -p workq -J pbv --mail-type=END,FAIL ./launch.pbvariants.smkj
```

If some jobs are failing, you will probably have to modify the content of the ```cluster.yaml``` to ask for more ressources.

If you believe the total job duration can exceed 4 days, use **-p unlimitq** or launch it on frontal, it spends most of it's time waiting and doesn't need much resources (see nohup bash command).

./launch.pbvariants.smkj -h will give a description of the options that are available. The option -u might be required if the run failed and was killed (out of memory, time limit exceeded) and has to be restarted.
