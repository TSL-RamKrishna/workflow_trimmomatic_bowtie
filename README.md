## Workflow for trimmomatic and bowtie
A snakemake workflow to run trimmomatic QC and bowtie mapping

## Requirements

1) snakemake (latest version) (source brew-default)
2) python v3.6+ (source python-3.6.1)
3) trimmomatic v0.36
4) bowtie2 v2.3.5


## How to run

Please source all the softwares first.

Run the following command in HPC
```
1) sbatch -o snakemake.log -J snakemake-workflow --wrap "snakemake -s analysis.smk analysis  --cluster-config cluster.json --cluster-config config.yaml  --cluster 'sbatch --mem {cluster.memory} --cpus {cluster.cpus}'  -p --jobs 4  --latency-wait 60 "  # job runs in parallel way in HPC
```

or on the local machine

```
snakemake -j 4 -s analysis.smk -p analysis
```

## How to run bowtie step and bam merge step indivitually

```
sbatch -o snakemake.log -J snakemake-workflow --wrap "snakemake -s analysis.smk run_bowtie  --cluster-config cluster.json --cluster-config config.yaml  --cluster 'sbatch --mem {cluster.memory} --cpus {cluster.cpus}'  -p --jobs 4  --latency-wait 60 " 

sbatch -o snakemake.log -J snakemake-workflow --wrap "snakemake -s analysis.smk run_merge_bam --cluster-config cluster.json --cluster-config config.yaml  --cluster 'sbatch --mem {cluster.memory} --cpus {cluster.cpus}'  -p --jobs 4  --latency-wait 60 " 
```
