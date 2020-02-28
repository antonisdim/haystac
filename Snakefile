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
        # todo find out why when they're are not comemnted out the rules from bowtie.smk doesn't work
        "example1/bowtie/example1.fasta",
        expand("example1/bowtie/example1.{n}.bt2", n=[1, 2, 3, 4]),
        expand("example1/bowtie/example1.rev.{n}.bt2", n=[1, 2]),
        "example1/bam_outputs/CSL.all.bam",
        "example1/bam_outputs/CSL.all_sorted.bam",
        "example1/bam_outputs/CSL.all_rmdup.bam",
        "example1/bam_outputs/CSL.all.fastq",
        "fastq/CSL.all.size"
