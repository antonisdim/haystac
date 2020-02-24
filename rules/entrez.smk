#!/usr/bin/env python
# -*- coding: utf-8 -*-

##### Target rules #####

rule entrez_nuccore_query:
    output:
        "entrez/{query}/{query}-nuccore.tsv"
    log:
        "entrez/{query}/{query}-nuccore.log"
    script:
        "../scripts/entrez_nuccore_query.py"

rule entrez_taxa_query:
    input:
        "entrez/{query}/{query}-nuccore.tsv"
    output:
        "entrez/{query}/{query}-taxa.tsv"
    log:
        "entrez/{query}/{query}-taxa.log"
    script:
        "../scripts/entrez_taxonomy_query.py"

checkpoint entrez_pick_sequences:
    input:
         "entrez/{query}/{query}-nuccore.tsv",
         "entrez/{query}/{query}-taxa.tsv"
    output:
         "entrez/{query}/{query}-selected-seqs.tsv"
    log:
         "entrez/{query}/{query}-selected-seqs.log"
    script:
        "../scripts/entrez_select_seqs.py"

rule entrez_download_sequence:
    log:
         "database/{orgname}/{accession}.log"
    output:
         "database/{orgname}/{accession}.fasta"
    script:
         "../scripts/entrez_download_sequence.py"

