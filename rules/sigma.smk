#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd


##### Target rules #####


rule index_database:
    input:
        "database/{orgname}/{accession}.fasta.gz"
    log:
        "database/{orgname}/{accession}_index.log"
    output:
        "database/{orgname}/{accession}_index.done"
    shell:
        "bowtie2-build --large-index {input} database/{wildcards.orgname}/{wildcards.accession} &> "
        "{log}; touch {output}"



MIN_FRAG_LEN = 15
MAX_FRAG_LEN = 150

rule alignments_per_taxon:
    input:
        fastq="{query}/fastq/{sample}_mapq.fastq.gz",
        bt2idx="database/{orgname}/{accession}.1.bt2l",
        readlen="{query}/fastq/{sample}_mapq.readlen"
    log:
        "{query}/sigma/{sample}/{orgname}/{accession}.log"
    output:
        "{query}/sigma/{sample}/{orgname}/{orgname}_{accession}.bam"
    params:
        # TODO use `lambda wildcards, input: ` to get the readlen file - Am I doing it wrong ?
        #  I need to read the file and then calculate the min score for the edit distance that's why it is like that
        min_score=lambda wildcards : (round(float(pd.read_csv(
            os.path.join(wildcards.query,'fastq',wildcards.sample + '_mapq.readlen'), header=None, squeeze=True) *
            float(config['mismatch_probability']))))*(-6) ,
        # min_score=-get_max_mismatch_count()*(-6),
        index="database/{orgname}/{accession}",
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
                  "-x {params.index} -U {input.fastq} | samtools sort -O bam -o {output} ) 2> {log}")
        elif config['sequencing'] == 'paired_end':
            #todo write the input files in the command correctly for paired end
            shell("(bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
                  "--score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {threads} "
                  "-I {params.minimum_fragment_length} -X {params.maximum_fragment_length} "
                  "-x {params.index} -U {input.fastq} | samtools sort -O bam -o {output} ) 2> {log}")

