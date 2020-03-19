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
        "{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam"
        # TODO index the BAM file here, rather than in the pysam code in parse_bams.py
    params:
        min_score=get_min_score,
        minimum_fragment_length=MIN_FRAG_LEN,
        maximum_fragment_length=MAX_FRAG_LEN
    threads:
        1
    run:
        # TODO split these into two separate rules -
        #  How can I call the right rule every time depending the sequencing strategy ?
        #  Shall I incorporate that in the output directories/names so that it knows what rule to use ?
        if config['sequencing'] == 'single_end':
            shell("(bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
                  "--score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {threads} "
                  "-x database/{wildcards.orgname}/{wildcards.accession} -U {input.fastq} | samtools sort -O bam -o {output} ) 2> {log}")
        elif config['sequencing'] == 'paired_end':
            #todo write the input files in the command correctly for paired end
            shell("(bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
                  "--score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {threads} "
                  "-I {params.minimum_fragment_length} -X {params.maximum_fragment_length} "
                  "-x database/{wildcards.orgname}/{wildcards.accession} -U {input.fastq} | samtools sort -O bam -o {output} ) 2> {log}")

