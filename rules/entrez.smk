#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

##### Target rules #####


checkpoint entrez_find_accessions:
    output:
        temp("{query}/entrez/{query}-accessions.tsv")
    log:
        temp("{query}/entrez/{query}-accessions.log")
    benchmark:
        repeat("benchmarks/entrez_find_accessions_{query}.benchmark.txt", 3)
    script:
        "../scripts/entrez_find_accessions.py"



rule entrez_nuccore_query:
    input:
        "{query}/entrez/{query}-accessions.tsv"
    output:
        temp("{query}/entrez/{query}_{chunk}-nuccore.tsv")
    log:
        temp("{query}/entrez/{query}_{chunk}-nuccore.log")
    benchmark:
        repeat("benchmarks/entrez_nuccore_query_{query}_{chunk}.benchmark.txt", 3)
    script:
        "../scripts/entrez_nuccore_query.py"



# noinspection PyUnresolvedReferences
def get_nuccore_chunks(wildcards):
    """
    Get all the accession chunks for the {query}-nuccore.tsv file.
    """

    chunk_size = 20  # todo Is this acceptable way to chunk the query based on chunk size ? Feel like I can do it better
    pick_accessions = checkpoints.entrez_find_accessions.get(query=wildcards.query)
    sequences = pd.read_csv(pick_accessions.output[0], sep='\t')

    if len(sequences) == 0:
        raise RuntimeError("The entrez find accessions file is empty.")

    if len(sequences) % chunk_size == 0:
        tot_chunks = len(sequences) / float(chunk_size)
    else:
        tot_chunks = (len(sequences) // float(chunk_size)) + 1

    inputs = []
    for chunk_num in range(int(tot_chunks)):
        inputs.append("{query}/entrez/{query}_{chunk}-nuccore.tsv".format(query=wildcards.query, chunk=chunk_num))

    return inputs



rule entrez_aggregate_nuccore:
    input:
        get_nuccore_chunks
    output:
        "{query}/entrez/{query}-nuccore.tsv"
    log:
        "{query}/entrez/{query}-nuccore.log"
    benchmark:
        repeat("benchmarks/entrez_aggregate_nuccore_{query}.benchmark.txt", 3)
    shell:
        "awk 'FNR>1 || NR==1' {input} 1> {output} 2> {log}"



rule entrez_taxa_query:
    input:
        "{query}/entrez/{query}-nuccore.tsv"
    output:
        "{query}/entrez/{query}-taxa.tsv"
    log:
        "{query}/entrez/{query}-taxa.log"
    benchmark:
        repeat("benchmarks/entrez_taxa_query_{query}.benchmark.txt", 3)
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
    benchmark:
        repeat("benchmarks/entrez_pick_sequences_{query}.benchmark.txt", 3)
    script:
        "../scripts/entrez_pick_sequences.py"



rule entrez_download_sequence:
    output:
        "database/{orgname}/{accession}.fasta.gz"
    log:
        "database/{orgname}/{accession}.log"
    benchmark:
        repeat("benchmarks/entrez_download_sequence_{orgname}_{accession}.benchmark.txt", 3)
    params:
        assembly=False
    script:
        "../scripts/entrez_download_sequence.py"



# noinspection PyUnresolvedReferences
def get_fasta_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep='\t')

    if len(sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

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
        "{query}/bowtie/{query}_entrez.fasta.gz"
    benchmark:
        repeat("benchmarks/entrez_multifasta_{query}.benchmark.txt", 3)
    shell:
        "cat {input} > {output}"
