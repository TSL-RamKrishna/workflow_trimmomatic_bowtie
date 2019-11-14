## Workflow for trimmomatic and bowtie
A snakemake workflow to run trimmomatic QC and bowtie mapping

## Requirements

1) snakemake (latest version) (source brew-default)
2) python v3.6+ (source python-3.6.1)
3) trimmomatic v0.36
4) bowtie2 v2.3.5

## How to run

Please source all the softwares first.

Different ways to run on the HPC (number 3 recommended)
```
1) snakemake -j 4 -s analysis.smk --cluster-config config.yaml --cluster-config cluster.json --cluster "sbatch analysis "  # runs in the headnode or interactive session, job gets killed if you close terminal/exit interactive.
2) sbatch -o logfile -J snakemake-workflow --wrap "snakemake -j 4 -s analysis.smk analysis -p "  # job runs in non parallel way, does not use HPC in optimized way
3) sbatch -o snakemake.log -J snakemake-workflow --wrap "snakemake -s analysis.smk analysis  --cluster-config cluster.json --cluster 'sbatch'  -p --jobs 4  --latency-wait 60 "  # job runs in parallel way in HPC
```

or on the local machine

```
snakemake -j 4 -s analysis.smk -p analysis
```
