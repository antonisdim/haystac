#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
from multiprocessing import cpu_count

##### Target rules #####

def get_fasta_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep='\t')

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = seq['GBSeq_organism'].replace(" ", "."), seq['GBSeq_accession-version']
        inputs.append('database/{orgname}/{accession}.fasta'.format(orgname=orgname, accession=accession))

    return inputs


rule bowtie_multifasta:
    input:
         get_fasta_sequences
    log:
         "{query}/bowtie/{query}.log"
    output:
         "{query}/bowtie/{query}.fasta"
    script:
          "../scripts/bowtie_multifasta.py"


rule bowtie_index:
    input:
         "{query}/bowtie/{query}.fasta"
    log:
         "{query}/bowtie/{query}.index.log"
    output:
         expand("{{query}}/bowtie/{{query}}.{n}.bt2", n=[1, 2, 3, 4]),
         expand("{{query}}/bowtie/{{query}}.rev.{n}.bt2", n=[1, 2])
    shell:
          "bowtie2-build {input} {wildcards.query}/bowtie/{wildcards.query} &> {log}"


rule bowtie_alignment:
    input:
        fastq=lambda wildcards: config['samples'][wildcards.sample],
        bt2idx="{query}/bowtie/{query}.1.bt2"
    log:
        "{query}/bam_outputs/{sample}.log"
    params:
        index="{query}/bowtie/{query}",
    output:
        "{query}/bam_outputs/{sample}_sorted.bam"
    threads:
        cpu_count()
    shell:
         "( bowtie2 -q --very-fast-local --threads {threads} -x {params.index} -U {input.fastq} "
         "| samtools sort -O bam -o {output} ) 2> {log}"


rule remove_duplicates:
    input:
        "{query}/bam_outputs/{sample}_sorted.bam"
    log:
        "{query}/bam_outputs/{sample}_rmdup.log"
    output:
        "{query}/bam_outputs/{sample}_rmdup.bam"
    params:
        folder=directory("{query}/bam_outputs/"),
        intermediate_file="{query}/bam_outputs/{sample}_sorted_rmdup.bam"
    shell:
        "dedup --merged --input {input} --output {params.folder}; "
        "mv {params.intermediate_file} {output} &> {log}"


rule extract_fastq:
    input:
        "{query}/bam_outputs/{sample}_rmdup.bam"
    log:
        "{query}/bam_outputs/{sample}_fastq.log"
    output:
        "{query}/bam_outputs/{sample}.fastq"
    params:
        config['alignment_qscore']
    shell:
        "samtools view -h -q {params} {input} | samtools fastq - > {output} &> {log}"


rule count_fastq_length:
    input:
         fastq=lambda wildcards: config['samples'][wildcards.sample]
    log:
         "fastq/{sample}.log"
    output:
         "fastq/{sample}.size"
    shell:
          "expr $(gunzip -c {input.fastq} | wc -l) / 4 1> {output} 2> {log}"

