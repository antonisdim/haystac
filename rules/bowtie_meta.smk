#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd

MIN_FRAG_LEN = 15
MAX_FRAG_LEN = 150

##### Target rules #####


rule index_database:
    input:
        "database/{orgname}/{accession}.fasta.gz"
    log:
        "database/{orgname}/{accession}_index.log"
    output:
        expand("database/{{orgname}}/{{accession}}.{n}.bt2l", n=[1, 2, 3, 4]),
        expand("database/{{orgname}}/{{accession}}.rev.{n}.bt2l", n=[1, 2])
    shell:
        "bowtie2-build --large-index {input} database/{wildcards.orgname}/{wildcards.accession} &> {log}"



def get_min_score(wildcards, input):
    return round(float(open(input.readlen).read()) * float(config['mismatch_probability'])) * -6  # TODO magic number -6



rule alignments_per_taxon:
    input:
        fastq="{query}/fastq/{sample}_mapq.fastq.gz",
        bt2idx="database/{orgname}/{accession}.1.bt2l",
        readlen="{query}/fastq/{sample}_mapq.readlen"
    log:
        "{query}/sigma/{sample}/{orgname}/{accession}.log"
    output:
        bam_file="{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam",
        bam_index="{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam.bai"
    params:
        min_score=get_min_score,
        minimum_fragment_length=MIN_FRAG_LEN,
        maximum_fragment_length=MAX_FRAG_LEN
    threads:
        1
    run:
        # TODO split these into two separate rules - don't know how
        if config['sequencing'] == 'single_end':
            shell("(bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
                  "--score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {threads} "
                  "-x database/{wildcards.orgname}/{wildcards.accession} -U {input.fastq} | "
                  "samtools sort -O bam -o {output.bam_file} ) 2> {log}; samtools index {output.bam_file}")

        elif config['sequencing'] == 'paired_end':
            #todo write the input files in the command correctly for paired end
            shell("(bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
                  "--score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {threads} "
                  "-I {params.minimum_fragment_length} -X {params.maximum_fragment_length} "
                  "-x database/{wildcards.orgname}/{wildcards.accession} -U {input.fastq} | "
                  "samtools sort -O bam -o {output} ) 2> {log}")


# noinspection PyUnresolvedReferences
def get_bamfile_paths(wildcards):
    """
    Get all the individual bam file paths for the taxa in our database.
    """
    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep='\t')

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = seq['species'].replace(" ", "_"), seq['GBSeq_accession-version']
        inputs.append('{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam'.
                      format(query=wildcards.query, sample=wildcards.sample, orgname=orgname, accession=accession))

    return inputs

# todo it is a bit of vile way of doing it but it is the only way of doing it just by calling rule_all,
#  happy to change it

rule all_alignments:
    input:
        get_bamfile_paths
    log:
        "{query}/sigma/{sample}_alignments.log"
    output:
        "{query}/sigma/{sample}_alignments.done"
    shell:
         "touch {output}"