#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from haystack.workflow.scripts.utilities import normalise_name


rule entrez_nuccore_query:
    output:
        config["db_output"] + "/entrez/entrez-nuccore.tsv",
    benchmark:
        repeat("benchmarks/entrez_nuccore_query.benchmark.txt", 1)
    message:
        "Fetching sequence metadata from the NCBI Nucleotide database for the query."
    script:
        "../scripts/entrez_nuccore_query.py"


rule entrez_taxa_query:
    input:
        config["db_output"] + "/entrez/entrez-nuccore.tsv",
    output:
        config["db_output"] + "/entrez/entrez-taxa.tsv",
    benchmark:
        repeat("benchmarks/entrez_taxa_query_entrez.benchmark.txt", 1)
    message:
        "Querying the NCBI Taxonomy database and fetching taxonomic metadata."
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_taxonomy_query.py"


def pick_after_refseq_prok(_):
    if config["refseq_rep"]:
        return config["db_output"] + "/entrez/genbank-genomes.tsv"
    else:
        return config["db_output"] + "/entrez/entrez-nuccore.tsv"


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
        "Selecting the longest sequence per taxon in the entrez query."
    script:
        "../scripts/entrez_pick_sequences.py"


rule entrez_download_sequence:
    output:
        config["cache"] + "/ncbi/{orgname}/{accession}.fasta.gz",
    benchmark:
        repeat("benchmarks/entrez_download_sequence_{orgname}_{accession}.benchmark.txt", 1)
    message:
        "Downloading accession {wildcards.accession} for taxon {wildcards.orgname}."
    wildcard_constraints:
        accession="[^-]+",
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_download_sequence.py"


def get_fasta_sequences(_):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    # noinspection PyUnresolvedReferences
    pick_sequences = checkpoints.entrez_pick_sequences.get()
    sequences = pd.read_csv(pick_sequences.output[0], sep="\t")

    if len(sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    inputs = []

    for key, seq in sequences.iterrows():
        orgname = normalise_name(seq["species"])
        accession = seq["AccessionVersion"]

        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz")

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
        "Concatenating all the fasta sequences for all the taxa of the entrez query."
    script:
        "../scripts/bowtie2_multifasta.py"
