#!/usr/bin/env python
# -*- coding: utf-8 -*-

configfile: "config.yaml"

##### Modules #####

include: "rules/entrez.smk"
include: "rules/bowtie.smk"
include: "rules/sigma.smk"

##### Target rules #####

rule all:
    input:
        "example1/entrez/example1-nuccore.tsv",
        "example1/entrez/example1-taxa.tsv",
        "example1/bowtie/example1.fasta",
        "example1/bowtie/example1.1.bt2",
        "example1/bam/CSL.all_sorted.bam",
        "example1/bam/CSL.all_sorted_rmdup.bam",
        "example1/bam/CSL.all.fastq",
        "fastq/CSL.all.size",
        "example1/bam/CSL.all.readlen",
        "database/Streptomyces.alboflavus/NZ_CP021748.1_index.done",
        "example1/sigma_outputs/CSL.all/Streptomyces.alboflavus/Streptomyces.alboflavus.bam"
