#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import pandas as pd

from haystack.workflow.scripts.utilities import normalise_name
from haystack.workflow.scripts.entrez_utils import get_accession_ftp_path


rule entrez_nuccore_query:
    output:
        config["db_output"] + "/entrez/entrez-nuccore.tsv",
    log:
        config["db_output"] + "/entrez/entrez-nuccore.tsv.log",
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
    log:
        config["db_output"] + "/entrez/entrez-taxa.log",
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


def get_rsync_url(wildcards):
    """Function to get NCBI urls for the database genomes"""

    try:
        url = get_accession_ftp_path(wildcards.accession, config)
        file_url = os.path.join(url, os.path.basename(url) + "_genomic.fna.gz")
        if file_url != "_genomic.fna.gz":
            return file_url
        else:
            return ""
    except RuntimeError:
        return ""
    except TypeError:
        # sometimes NCBI returns a None type url, but the URL does exist if I do it independently
        get_rsync_url(wildcards)


rule entrez_download_sequence:
    output:
        config["cache"] + "/{orgname}/{accession}.fasta.gz",
    log:
        config["cache"] + "/{orgname}/{accession}.fasta.gz.log",
    benchmark:
        repeat("benchmarks/entrez_download_sequence_{orgname}_{accession}.benchmark.txt", 1)
    params:
        url=get_rsync_url,
    message:
        "Downloading accession {wildcards.accession} for taxon {wildcards.orgname}."
    wildcard_constraints:
        accession="[^-]+",
    shell:
        "(if [[ -n '{params.url}' && ! `{config[mtDNA]}` ]]; "
        "   then (wget -q -O - '{params.url}' | gunzip ); "
        "   else (python {config[workflow_dir]}/scripts/entrez_download_sequence.py"
        "           --database nuccore "
        "           --accession {wildcards.accession} ); "
        "fi) | bgzip -f 1> {output} 2> {log}"


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

        inputs.append(config["cache"] + "/{orgname}/{accession}.fasta.gz".format(orgname=orgname, accession=accession))

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
