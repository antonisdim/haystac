#!/usr/bin/env python
# -*- coding: utf-8 -*-

configfile: "config.yaml"

##### Modules #####

include: "rules/entrez.smk"

##### Target rules #####

rule all:
    input:
        "entrez/example1/example1-nuccore.tsv",
        "entrez/example1/example1-taxa.tsv",
        "entrez/example1/example1-selected-seqs.tsv",
        "pipeline_files/example1/refseq/",
        "pipeline_files/example1/reference_genome/example1.fasta"