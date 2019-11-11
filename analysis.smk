#!/usr/bin/env python3

configfile: "config.yaml"

#print(config['Samples']['sampleA']['R1'])
#print(config['Samples']['sampleB']['R1'])

reference=config['Reference']
analysis_samples = config['Samples'].keys()
print('samples are ', list(analysis_samples))
# for sample in analysis_samples:
#     print(sample)



rule trimmomatic:
    input:
        R1=lambda wildcards: config['Samples'][wildcards.sample]['R1'],
        R2=lambda wildcards: config['Samples'][wildcards.sample]['R2']
        # R1=config['Samples'][sample]['R1'],
        # R2=config['Samples'][sample]['R2']
    output:
        R1Paired="results/trimmomatic/{sample}/R1_paired.fastq",
        R1Unpaired="results/trimmomatic/{sample}/R1_unpaired.fastq",
        R2Paired="results/trimmomatic/{sample}/R2_paired.fastq",
        R2Unpaired="results/trimmomatic/{sample}/R2_unpaired.fastq"
    message:
        'Running Quality Control using Trimmomatic '
    shell:
        "source trimmomatic-0.36; trimmomatic PE -threads 4  -phred64 -trimlog logs/{wildcards.sample}/trimmomatic.log -validatePairs {input.R1} {input.R2} {output} SLIDINGWINDOW:4:20 MINLEN:30"

rule bowtie_indexing:
    input: config['Reference']
    output: config["Reference"] + ".1.bt2", config["Reference"] + ".2.bt2", config["Reference"] + ".3.bt2", config["Reference"] + ".4.bt2", config["Reference"] + ".rev.1.bt2", config["Reference"] + ".rev.2.bt2"
    shell:
        "bowtie2-build -f {input} {input}"

rule bowtie_mapping:
    input:
        R1="results/trimmomatic/{sample}/R1_paired.fastq",
        R2="results/trimmomatic/{sample}/R2_paired.fastq"
        index=[config["Reference"] + ".1.bt2",
	config["Reference"] + ".2.bt2",
	config["Reference"] + ".3.bt2",
	config["Reference"] + ".4.bt2",
	config["Reference"] + ".rev.1.bt2",
	config["Reference"] + ".rev.2.bt2"]
    output:
        temp("results/bowtie_mapping/{sample}/alignment.sam")

    shell:
        "source bowtie2-2.3.5; bowtie2 -x {input.index} --maxins 1000  -1 {input.R1} -2 {input.R2} -S {output}"

rule samtobam:
	input: "results/bowtie_mapping/{sample}/alignment.sam"
	output: temp("results/bowtie_mapping/{sample}/alignment.bam")
	shell: "source samtools-1.9; samtools view -b -T {reference} -o {output} {input}"

rule sort_bam:
    input: "results/bowtie_mapping/{sample}/alignment.bam"
    output: protected("results/bowtie_mapping/{sample}/alignment_sorted.bam")
    shell: "source samtools-1.9; samtools sort -l 5 -o {output} {input}"   # -l option is compress level.

rule analysis:
    input:
        expand("results/bowtie_mapping/{sample}/alignment_sorted.bam",
                sample=analysis_samples)
