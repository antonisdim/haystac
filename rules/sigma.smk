#!/usr/bin/env python
# -*- coding: utf-8 -*-

##### Target rules #####

import pandas as pd
import os

rule index_database:
    input:
        "database/{orgname}/{accession}.fasta"
    log:
        "database/{orgname}/{accession}_index.log"
    output:
        # todo doing it in a dirty way and not as in the commented one,
        #  cause otherwise I've got issues of wildcards not matching in this and the next rule.
        #  If this alright I'll delete the commented lines.
        # expand("database/{{orgname}}/{{orgname}}.{n}.bt2", n=[1, 2, 3, 4]),
        # expand("database/{{orgname}}/{{orgname}}.rev.{n}.bt2", n=[1, 2])
        "database/{orgname}/{accession}_index.done"
    shell:
        "bowtie2-build {input} database/{wildcards.orgname}/{wildcards.orgname} &> {log}; touch {output}"


rule alignments_per_taxon:
    input:
        fastq="{query}/bam/{sample}.fastq",
        bt2idx="database/{orgname}/{orgname}.1.bt2",
        readlen="{query}/bam/{sample}.readlen"
    log:
        "{query}/sigma_outputs/{sample}/{orgname}/{orgname}.log"
    output:
        "{query}/sigma_outputs/{sample}/{orgname}/{orgname}.bam",
    params:
        min_score=lambda wildcards : (round(float(pd.read_csv(
            os.path.join(wildcards.query,'bam',wildcards.sample + '.readlen'), header=None, squeeze=True) *
            float(config['mismatch_probability']))))*(-6) ,
        # min_score=-get_max_mismatch_count()*(-6),
        index="database/{orgname}/{orgname}",
        bowtie_threads_number=1,
        minimum_fragment_length=15,
        maximum_fragment_length=150
    run:
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

