#!/usr/bin/env python
# -*- coding: utf-8 -*-

configfile: "config.yaml"

##### Modules #####

include: "rules/entrez.smk"
include: "rules/bowtie.smk"
include: "rules/bowtie_meta.smk"
include: "rules/metagenomics.smk"
include: "rules/entrez_alt.smk"

##### Wildcards #####

wildcard_constraints:
    query="[\w]+",
    sample="[\w]+",
    # orgname="[\w]+", #todo these wildcards need to have numbers, letters and symbols like '_, .'
    # accession="[\w]+",
    chunk='\d+'

##### Target rules #####

from datetime import datetime
startTime = datetime.now()
rule all:
    input:
        "yersinia_test/entrez_alt/sizes.txt",
        # "yersinia_test/entrez/yersinia_test-nuccore.tsv",
        # "yersinia_test/bowtie/yersinia_test.fasta.gz",  # test entrez.smk
        # "yersinia_test/fastq/RISE00_mapq.readlen",      # test bowtie.smk
        # "yersinia_test/sigma/RISE00_alignments.done",   # test bowtie_meta.smk
        # "yersinia_test/probabilities/RISE00/RISE00_posterior_probabilities.csv", #test metagenomics.smk - probabilities
        # "yersinia_test/probabilities/RISE00/RISE00_posterior_abundance.tsv" #test metagenomics.smk - abundances

print(datetime.now() - startTime)


