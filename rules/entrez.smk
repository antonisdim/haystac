#!/usr/bin/env python
# -*- coding: utf-8 -*-

##### Target rules #####

rule entrez_nuccore_query:
    output:
        "{query}/entrez/{query}-nuccore.tsv"
    log:
        "{query}/entrez/{query}-nuccore.log"
    script:
        "../scripts/entrez_nuccore_query.py"


rule entrez_taxa_query:
    input:
        "{query}/entrez/{query}-nuccore.tsv"
    output:
        "{query}/entrez/{query}-taxa.tsv"
    log:
        "{query}/entrez/{query}-taxa.log"
    script:
        "../scripts/entrez_taxonomy_query.py"

# todo it doesn't work as a checkpoint, even when I call a top level rule like the bowtie2_multifasta one
checkpoint entrez_pick_sequences:
    input:
         "{query}/entrez/{query}-nuccore.tsv",
         "{query}/entrez/{query}-taxa.tsv"
    output:
         "{query}/entrez/{query}-selected-seqs.tsv"
    log:
         "{query}/entrez/{query}-selected-seqs.log"
    script:
        "../scripts/entrez_pick_sequences.py"

# todo doesn't work
rule entrez_download_sequence:
    output:
        "database/{orgname}/{accession}.fasta"
    log:
        "database/{orgname}/{accession}.log"
    script:
         "../scripts/entrez_download_sequence.py"

