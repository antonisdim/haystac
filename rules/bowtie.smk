#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import os

##### Target rules #####

def get_fasta_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep='\t')

    inputs = []

    for key, seq in sequences.iterrows():
        inputs.append('database/{orgname}/{accession}.fasta'.format(orgname=seq['TSeq_orgname'].replace(" ", "."),
                                                                    accession=seq['TSeq_accver']))

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

# todo don't know how to run rules with multiple outputs
rule bowtie_index:
    input:
         "{query}/bowtie/{query}.fasta"
    log:
         "{query}/bowtie/{query}.bt2.log"
    output:
         expand("{{query}}/bowtie/{{query}}.{n}.bt2", n=[1, 2, 3, 4]),
         expand("{{query}}/bowtie/{{query}}.rev.{n}.bt2", n=[1, 2])
    shell:
          "bowtie2-build {input} {wildcards.query}/bowtie/{wildcards.query} &> {log}"


def trim_ext(fullpath, n=1):
    return ('.').join(fullpath.split('.')[:-n])

def trim_path_ext(fullpath, N=1):
    return trim_ext(os.path.basename(fullpath), n=N)


rule bowtie_alignment:
    input:
         # todo if i use just {sample} as a wildcard I get and error of it being unresovled, why ? Why does the query wildcard just work ?
         lambda wildcards: config['samples'][wildcards.sample]
    log:
       # todo is there a way to trim the sample name from the sample wildcord on the fly ?
         "{query}/bam_outputs/{trim_sample}.log"
    output:
        "{query}/bam_outputs/{trim_sample}.bam"
    shell:
          "bowtie2 -q --very-fast-local -p 20 -x {wildcards.query}/bowtie/{wildcards.query} -U {input} " \
          "| samtools view -Shu > {output}"



rule count_fastq_length:
    input:
         lambda wildcards: config['samples'][wildcards.sample]
    log:
         "fastq/{sample}.log"
    output:
         "fastq/{sample}.size"
    shell:
          "expr $(gunzip -c {input} | wc -l) / 4 1> {output} 2> {log}"

