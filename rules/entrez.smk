#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import sys
import pandas as pd


##### Target rules #####


from scripts.entrez_nuccore_query import CHUNK_SIZE
from scripts.rip_utilities import normalise_name


checkpoint entrez_find_accessions:
    output:
        temp(config["db_output"] + "/entrez/entrez-accessions.tsv"),
    log:
        temp(config["db_output"] + "/entrez/entrez-accessions.log"),
    benchmark:
        repeat("benchmarks/entrez_find_accessions.benchmark.txt", 1)
    message:
        "Finding all the accessions, whose metadata are going to be fetched, for the entrez query. "
        "The temporary output can be found in {output} and the its log file in {log}."
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_find_accessions.py"


rule entrez_nuccore_query:
    input:
        config["db_output"] + "/entrez/entrez-accessions.tsv",
    output:
        temp(config["db_output"] + "/entrez/entrez_{chunk}-nuccore.tsv"),
    log:
        config["db_output"] + "/entrez/entrez_{chunk}-nuccore.log",
    benchmark:
        repeat("benchmarks/entrez_nuccore_query_entrez_{chunk}.benchmark.txt", 1)
    message:
        "Fetching sequence metadata from the NCBI Nucleotide database "
        "for the accessions in chunk {wildcards.chunk} for the NCBI entrez query. "
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

    pick_accessions = checkpoints.entrez_find_accessions.get()
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
            config["db_output"]
            + "/entrez/entrez_{chunk}-nuccore.tsv".format(chunk=chunk_num)
        )

    return inputs


rule entrez_aggregate_nuccore:
    input:
        get_nuccore_chunks,
    output:
        config["db_output"] + "/entrez/entrez-nuccore.tsv",
    benchmark:
        repeat("benchmarks/entrez_aggregate_nuccore_entrez.benchmark.txt", 1)
    message:
        "Concatenating all the temporary output files containing accession metadata fetched from the NCBI Nucleotide "
        "database for the entrez query. The permanent output can be found in {output}."
    shell:
        "awk 'FNR>1 || NR==1' {input} 1> {output}"


rule entrez_taxa_query:
    input:
        config["db_output"] + "/entrez/entrez-nuccore.tsv",
    output:
        config["db_output"] + "/entrez/entrez-taxa.tsv",
    log:
        config["db_output"] + "/entrez/entrez-taxa.log",
    benchmark:
        repeat("benchmarks/entrez_taxa_query_entrez.benchmark.txt", 1)
    message:
        "Querying the NCBI Taxonomy database to fetch taxonomic metadata, for all the different "
        "taxa the NCBI Nucleotide returned accessions belong to, for the entrez query. "
        "The output table can be found in {output} and its log file in {log}."
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_taxonomy_query.py"


def pick_after_refseq_prok(wildcards):

    if config["refseq_rep"]:
        return config["db_output"] + "/entrez/genbank-genomes.tsv"
    else:
        return config["db_output"] + "/entrez/entrez_1-nuccore.tsv"


checkpoint entrez_pick_sequences:
    input:
        nuccore=config["db_output"] + "/entrez/entrez-nuccore.tsv",
        taxonomy=config["db_output"] + "/entrez/entrez-taxa.tsv",
        priority=pick_after_refseq_prok,
    output:
        config["db_output"] + "/entrez/entrez-selected-seqs.tsv",
    log:
        config["db_output"] + "/entrez/entrez-selected-seqs.log",
    benchmark:
        repeat("benchmarks/entrez_pick_sequences_entrez.benchmark.txt", 1)
    message:
        "Selecting the longest sequence per taxon for entrez query. "
        "Input tables with accession and taxonomic metadata can be found in "
        "{input.nuccore} and {input.taxonomy} respectively."
        "The output table can be found in {output} and its log file in {log}. "
    script:
        "../scripts/entrez_pick_sequences.py"


rule entrez_download_sequence:
    output:
        config["genome_cache_folder"] + "/{orgname}/{accession}.fasta.gz",
    log:
        config["genome_cache_folder"] + "/{orgname}/{accession}.log",
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
        orgname = normalise_name(seq["species"])
        accession = seq["GBSeq_accession-version"]

        inputs.append(
            config["genome_cache_folder"]
            + "/{orgname}/{accession}.fasta.gz".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


rule entrez_multifasta:
    input:
        get_fasta_sequences,
    log:
        config["db_output"] + "/bowtie/entrez_query.log",
    output:
        config["db_output"] + "/bowtie/entrez_query.fasta.gz",
    benchmark:
        repeat("benchmarks/entrez_multifasta_entrez_query.benchmark.txt", 1)
    message:
        "Concatenating all the fasta sequences for all the taxa in {output} for the entrez query, and its "
        "log file can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"
