## Workflow for trimmomatic and bowtie
A snakemake workflow to run trimmomatic QC and bowtie mapping

## Requirements

1) python v3.7+
2) snakemake (latest version)
3) trimmomatic v0.36
4) bowtie2 v2.3.5

## How to run

Run on the HPC
```
1) snakemake -j 4 -s analysis.smk --cluster-config config.yaml --cluster-config cluster.json --cluster "sbatch analysis "
2) source snakemake-5.7.4; sabtch -o logfile -J snakemake-workflow --wrap "snakemake -j 4 -s analysis.smk analysis -p "
```

or on the local machine

```
snakemake -j 4 -s analysis.smk -p analysis
```
