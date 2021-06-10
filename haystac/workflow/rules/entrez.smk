#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from haystac.workflow.scripts.utilities import get_total_paths


rule entrez_nuccore_query:
    output:
        config["db_output"] + "/entrez/entrez-nuccore.tsv",
    message:
        "Fetching sequence metadata from the NCBI Nucleotide database for the query."
    script:
        "../scripts/entrez_nuccore_query.py"


rule entrez_taxa_query:
    input:
        config["db_output"] + "/entrez/entrez-nuccore.tsv",
    output:
        config["db_output"] + "/entrez/entrez-taxa.tsv",
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
    message:
        "Selecting the longest sequence per taxon in the entrez query."
    script:
        "../scripts/entrez_pick_sequences.py"


rule entrez_download_sequence:
    output:
        config["cache"] + "/ncbi/{orgname}/{accession}.fasta.gz",
    message:
        "Downloading accession {wildcards.accession} for taxon {wildcards.orgname}."
    wildcard_constraints:
        accession="[^-]+",
    resources:
        entrez_api=1,
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/entrez_download_sequence.py"


def get_all_accessions(_):
    """Get fasta paths in our db"""
    return [
        config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz"
        for orgname, accession in get_total_paths(checkpoints, config)
    ]


checkpoint entrez_db_list:
    input:
        get_all_accessions,
    log:
        config["db_output"] + "/db_taxa_accessions.log",
    output:
        config["db_output"] + "/db_taxa_accessions.tsv",
    message:
        "Aggregating all the species/accession pairs that exist in the database."
    script:
        "../scripts/entrez_db_list.py"
