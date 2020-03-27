#!/usr/bin/env python
# -*- coding: utf-8 -*-

configfile: "config.yaml"

##### Modules #####

include: "rules/entrez.smk"
include: "rules/bowtie.smk"
include: "rules/bowtie_meta.smk"
include: "rules/metagenomics.smk"

##### Target rules #####

rule all:
    input:
        "yersinia_test/entrez/yersinia_test-nuccore.tsv",
        "yersinia_test/bowtie/yersinia_test.fasta.gz",  # test entrez.smk
        "yersinia_test/fastq/RISE00_mapq.readlen",      # test bowtie.smk
        "yersinia_test/sigma/RISE00_alignments.done",   # test bowtie_meta.smk
        "yersinia_test/probabilities/RISE00/RISE00_posterior_abundance.tsv"




