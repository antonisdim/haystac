#!/usr/bin/env python
# -*- coding: utf-8 -*-

configfile: "config.yaml"

##### Modules #####

include: "rules/entrez.smk"
include: "rules/bowtie.smk"

##### Target rules #####

rule all:
    input:
        # "example1/entrez/example1-nuccore.tsv",
        # "example1/entrez/example1-taxa.tsv",
        # "example1/entrez/example1-selected-seqs.tsv",
        # todo if the rules from entrez.smk are not commented out the rules from the bowtie.smk don't get executed
        "example1/bowtie/example1.fasta",
        # "example1/bowtie/example1.1.bt2"
        # "example1/bam_outputs/CSL.all.bam",
        # "example1/bam_outputs/CSL.all_sorted.bam",
        # "example1/bam_outputs/CSL.all_rmdup.bam",
        # "example1/bam_outputs/CSL.all.fastq",
        # "fastq/CSL.all.size"
