#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

sys.path.append(os.getcwd())
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))


##### Target rules #####

WITH_REFSEQ_REP = config["WITH_REFSEQ_REP"]

from scripts.entrez_nuccore_query import CHUNK_SIZE


checkpoint entrez_find_accessions:
    output:
        temp("{query}/entrez/{query}-accessions.tsv"),
    log:
        temp("{query}/entrez/{query}-accessions.log"),
    benchmark:
        repeat("benchmarks/entrez_find_accessions_{query}.benchmark.txt", 1)
    message:
        "Finding all the accessions, whose metadata are going to be fetched, for query {wildcards.query}. "
        "The temporary output can be found in {output} and the its log file in {log}."
    script:
        "../scripts/entrez_find_accessions.py"


rule entrez_nuccore_query:
    input:
        "{query}/entrez/{query}-accessions.tsv",
    output:
        temp("{query}/entrez/{query}_{chunk}-nuccore.tsv"),
    log:
        temp("{query}/entrez/{query}_{chunk}-nuccore.log"),  # TODO why are the logs temporary?
    benchmark:
        repeat("benchmarks/entrez_nuccore_query_{query}_{chunk}.benchmark.txt", 1)
    message:
        # TODO it looks weird that the scheduler picks the chunks at random
        #   use the `priority` attribute to force them to download in order
        #   see https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#priorities
        "Fetching sequence metadata from the NCBI Nucleotide database "
        "for the accessions in chunk {wildcards.chunk} for query {wildcards.query}. "
        "The temporary output can be found in {output} and its log file in {log}."
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_nuccore_query.py"


# noinspection PyUnresolvedReferences
def get_nuccore_chunks(wildcards):
    """
    Get all the accession chunks for the {query}-nuccore.tsv file.
    """

    pick_accessions = checkpoints.entrez_find_accessions.get(query=wildcards.query)
    sequences = pd.read_csv(pick_accessions.output[0], sep="\t")

    if len(sequences) == 0:
        raise RuntimeError("The entrez find accessions file is empty.")

    if len(sequences) % CHUNK_SIZE == 0:
        tot_chunks = len(sequences) / float(CHUNK_SIZE)
    else:
        tot_chunks = (len(sequences) // float(CHUNK_SIZE)) + 1

    inputs = []
    for chunk_num in range(int(tot_chunks)):
        inputs.append(
            "{query}/entrez/{query}_{chunk}-nuccore.tsv".format(
                query=wildcards.query, chunk=chunk_num
            )
        )

    return inputs


rule entrez_aggregate_nuccore:
    input:
        get_nuccore_chunks,
    output:
        "{query}/entrez/{query}-nuccore.tsv",
    log:
        "{query}/entrez/{query}-nuccore.log",
    benchmark:
        repeat("benchmarks/entrez_aggregate_nuccore_{query}.benchmark.txt", 1)
    message:
        "Concatenating all the temporary output files containing accession metadata fetched from the NCBI Nucleotide "
        "database for query {wildcards.query}. The permanent output can be found in {output} and its log file in {log}."
    shell:
        "awk 'FNR>1 || NR==1' {input} 1> {output} 2> {log}"


rule entrez_taxa_query:
    input:
        "{query}/entrez/{query}-nuccore.tsv",
    output:
        "{query}/entrez/{query}-taxa.tsv",
    log:
        "{query}/entrez/{query}-taxa.log",
    benchmark:
        repeat("benchmarks/entrez_taxa_query_{query}.benchmark.txt", 1)
    message:
        "Querying the NCBI Taxonomy database to fetch taxonomic metadata, for all the different "
        "taxa the NCBI Nucleotide returned accessions belong to, for query {wildcards.query}. "
        "The output table can be found in {output} and its log file in {log}."
    script:
        "../scripts/entrez_taxonomy_query.py"


def pick_after_refseq_prok(wildcards):

    if WITH_REFSEQ_REP:
        return "{query}/entrez/{query}-genbank-genomes.tsv".format(
            query=wildcards.query
        )
    else:
        return "{query}/entrez/{query}_1-nuccore.tsv"


checkpoint entrez_pick_sequences:
    input:
        nuccore="{query}/entrez/{query}-nuccore.tsv",
        taxonomy="{query}/entrez/{query}-taxa.tsv",
        priority=pick_after_refseq_prok,
    output:
        "{query}/entrez/{query}-selected-seqs.tsv",
    log:
        "{query}/entrez/{query}-selected-seqs.log",
    benchmark:
        repeat("benchmarks/entrez_pick_sequences_{query}.benchmark.txt", 1)
    message:
        "Selecting the longest sequence per taxon for query {wildcards.query}. "
        "Input tables with accession and taxonomic metadata can be found in "
        "{input.nuccore} and {input.taxonomy} respectively."
        "The output table can be found in {output} and its log file in {log}. "
    script:
        "../scripts/entrez_pick_sequences.py"


rule entrez_download_sequence:
    output:
        "database/{orgname}/{accession}.fasta.gz",
    log:
        "database/{orgname}/{accession}.log",
    benchmark:
        repeat("benchmarks/entrez_download_sequence_{orgname}_{accession}.benchmark.txt", 1)
    params:
        assembly=False,
    wildcard_constraints:
        accession="\w+\.\d+", # TODO refactor this so we're not reliant on the style of the accession (low priority)
    message:
        "Downloading accession {wildcards.accession} for taxon {wildcards.orgname}. "
        "The downloaded fasta sequence can be found in {output} and its log file in {log}."
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_download_sequence.py"


# noinspection PyUnresolvedReferences
def get_fasta_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
    sequences = pd.read_csv(pick_sequences.output[0], sep="\t")

    if len(sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
            seq["GBSeq_accession-version"],
        )
        inputs.append(
            "database/{orgname}/{accession}.fasta.gz".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


rule entrez_multifasta:
    input:
        get_fasta_sequences,
    log:
        "{query}/bowtie/{query}.log",
    output:
        "{query}/bowtie/{query}_entrez.fasta.gz",
    benchmark:
        repeat("benchmarks/entrez_multifasta_{query}.benchmark.txt", 1)
    message:
        "Concatenating all the fasta sequences for all the taxa in {output} for query {wildcards.query}, and its "
        "log file can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"
