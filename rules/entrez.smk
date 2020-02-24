#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

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


def get_fasta_sequences(wildcards):

    seqs = checkpoints.entrez_pick_sequences.get(query=wildcards.query)

    sequences = pd.read_csv(seqs.output[0], sep='\t')

    inputs = []

    for key, seq in sequences.iterrows():
        orgname = seq['TSeq_orgname'].replace(" ", ".")
        accession = seq['TSeq_accver']

        inputs.append('database/{orgname}/{accession}.fasta'.format(orgname=orgname, accession=accession))

    return inputs


rule make_bowtie_general_profile:
    input:
         get_fasta_sequences
    log:
         "bowtie/{query}/{query}.log"
    output:
         "bowtie/{query}/{query}.fasta"
    script:
          "../scripts/entrez_bowtie_multifasta.py"


