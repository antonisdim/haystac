#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

##### Target rules #####

# TODO this is so slow... what if we split this into three rules?
#      1. run query - checkpoint rule, w/ temporary() output containing the
#      2. fetch chunk - can be run in parallel, with delayed start to avoid 3 query per sec limit
#                       e.g. if we pass in a chunk_num param, and use `sleep(chunk_num // 3)`
#      3. join all chunks

checkpoint entrez_find_accesions:
    output:
        temp("{query}/entrez/{query}-accessions.tsv")
    log:
        temp("{query}/entrez/{query}-accessions.log")
    script:
        "../scripts/entrez_find_accessions.py"



rule entrez_nuccore_query:
    input:
        "{query}/entrez/{query}-accessions.tsv"
    output:
        temp("{query}/entrez/{chunk}-nuccore.tsv")
    log:
        temp("{query}/entrez/{chunk}-nuccore.log")
    script:
        "../scripts/entrez_nuccore_query.py"



# noinspection PyUnresolvedReferences
def get_nuccore_chunks(wildcards):
    """
    Get all the accession chunks for the {query}-nuccore.tsv file.
    """

    CHUNKS = 50 #todo maybe this should be in the cofig ?
    pick_accessions = checkpoints.entrez_find_accesions.get(query=wildcards.query)
    sequences = pd.read_csv(pick_accessions.output[0], sep='\t')

    inputs = []

    for chunk_num in range(CHUNKS):
        inputs.append("{query}/entrez/{chunk}-nuccore.tsv".format(query=wildcards.query, chunk=chunk_num))

    return inputs



rule entrez_aggregate_nuccore:
    input:
        get_nuccore_chunks
    output:
        "{query}/entrez/{query}-nuccore.tsv"
    log:
        "{query}/entrez/{query}-nuccore.log"
    shell:
        'cat {input} 1> {output} 2> {log}'



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
    output:
        "database/{orgname}/{accession}.fasta.gz"
    log:
        "database/{orgname}/{accession}.log"
    script:
         "../scripts/entrez_download_sequence.py"



# noinspection PyUnresolvedReferences
def get_fasta_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep='\t')

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = seq['species'].replace(" ", "_"), seq['GBSeq_accession-version']
        inputs.append('database/{orgname}/{accession}.fasta.gz'.format(orgname=orgname, accession=accession))

    return inputs



rule entrez_multifasta:
    input:
         get_fasta_sequences
    log:
         "{query}/bowtie/{query}.log"
    output:
         "{query}/bowtie/{query}.fasta.gz"
    shell:
         "cat {input} > {output}"
