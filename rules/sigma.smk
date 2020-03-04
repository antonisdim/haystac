#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd


##### Target rules #####


rule index_database:
    input:
        "database/{orgname}/{accession}.fasta"
    log:
        "database/{orgname}/{accession}_index.log"
    output:
        # TODO go back to the old way of doing this, and use the {orgname}/{accession} wildcards, not {orgname}}/{orgname}
        #      this is important because if a second search is run, then the longest sequence may have changed, in which
        #      case the index needs to be rebuilt
        # expand("database/{{orgname}}/{{orgname}}.{n}.bt2l", n=[1, 2, 3, 4]),
        # expand("database/{{orgname}}/{{orgname}}.rev.{n}.bt2l", n=[1, 2])
        "database/{orgname}/{accession}_index.done"
    shell:
         # TODO did you forget to tell bowtie2 to use a --large-index?
        "bowtie2-build {input} database/{wildcards.orgname}/{wildcards.orgname} &> {log}; touch {output}"


rule alignments_per_taxon:
    input:
        fastq="{query}/fastq/{sample}_mapq.fastq",
        bt2idx="database/{orgname}/{orgname}.1.bt2l",
        readlen="{query}/bam/{sample}.readlen"
    log:
        "{query}/sigma/{sample}/{orgname}/{orgname}.log"
    output:
        "{query}/sigma/{sample}/{orgname}/{orgname}.bam",
    params:
        # TODO use `lambda wildcards, input: ` to get the readlen file
        min_score=lambda wildcards, input : (round(float(pd.read_csv(
            os.path.join(wildcards.query,'bam',wildcards.sample + '.readlen'), header=None, squeeze=True) *
            float(config['mismatch_probability']))))*(-6) ,
        # min_score=-get_max_mismatch_count()*(-6),
        index="database/{orgname}/{orgname}",
        bowtie_threads_number=1,     # TODO use `threads` not `params` for this
        minimum_fragment_length=15,  # TODO no magic numbers, use static constants
        maximum_fragment_length=150
    run:
        # TODO split these into two separate rules
        if config['sequencing'] == 'single_end':
            shell("(bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
                  "--score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {params.bowtie_threads_number} "
                  "-x {params.index} -U {input.fastq} | samtools sort -O bam -o {output} ) 2> {log}")
        elif config['sequencing'] == 'paired_end':
            #todo write the input files in the command correctly for paired end
            shell("(bowtie2 --time --no-unal --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 "
                  "--score-min L,{params.min_score},0.0 --gbar 1000 -q --threads {params.bowtie_threads_number} "
                  "-I {params.minimum_fragment_length} -X {params.maximum_fragment_length} "
                  "-x {params.index} -U {input.fastq} | samtools sort -O bam -o {output} ) 2> {log}")

