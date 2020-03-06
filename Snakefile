#!/usr/bin/env python
# -*- coding: utf-8 -*-

configfile: "config.yaml"

##### Modules #####

include: "rules/entrez.smk"
include: "rules/bowtie.smk"
include: "rules/sigma.smk"  # TODO rename file now that we're not using sigma

##### Target rules #####

rule all:
    input:
        "yersinia_test/entrez/yersinia_test-nuccore.tsv",
        "yersinia_test/entrez/yersinia_test-taxa.tsv",
        "yersinia_test/bowtie/yersinia_test.fasta",
        "yersinia_test/bowtie/yersinia_test.1.bt2l",
        "yersinia_test/bam/RISE00_sorted.bam",
        "yersinia_test/bam/RISE00_sorted_rmdup.bam",
        # "yersinia_test/bam/RISE00_mapq.fastq",
        # "fastq/RISE00.size",
        # "yersinia_test/bam/RISE00.readlen",
        # "database/Streptomyces.alboflavus/NZ_CP021748.1_index.done",
        # "example1/sigma_outputs/CSL.all/Streptomyces.alboflavus/Streptomyces.alboflavus.bam"
