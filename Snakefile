#!/usr/bin/env python
# -*- coding: utf-8 -*-

configfile: "config.yaml"

##### Modules #####

include: "rules/entrez.smk"
include: "rules/bowtie.smk"

##### Target rules #####

rule all:
    input:
        # todo find out why when they're are not comemnted out the rule from bowtie.smk doesn't work
        # "example1/entrez/example1-nuccore.tsv",
        # "example1/entrez/example1-taxa.tsv",
        # "example1/entrez/example1-selected-seqs.tsv",
        "example1/bowtie/example1.fasta",
         "example1/bam_outputs/CSL.all.bam"
