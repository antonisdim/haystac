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

rule entrez_download_sequence:
    log:
         "sequences/{orgname}/{accession}.log"
    output:
         "sequences/{orgname}/{accession}.fasta"
    script:
         "../scripts/entrez_download_sequence.py"

