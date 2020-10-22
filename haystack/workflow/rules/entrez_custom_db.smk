#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from haystack.workflow.scripts.utilities import normalise_name, check_unique_taxa_in_custom_input


rule entrez_custom_sequences:
    input:
        config["sequences"],
    log:
        config["cache"] + "/ncbi/{orgname}/custom_seq-{accession}.log",
    output:
        config["cache"] + "/ncbi/{orgname}/custom_seq-{accession}.fasta.gz",
    message:
        "Adding the user provided fasta sequence {wildcards.accession} for taxon {wildcards.orgname} to the database."
    script:
        "../scripts/entrez_custom_sequences.py"


def get_paths_for_custom_seqs():
    if config["sequences"] == "":
        return ""

    custom_fasta_paths = pd.read_csv(
        config["sequences"], sep="\t", header=None, names=["species", "accession", "path"],
    )

    if len(custom_fasta_paths) == 0:
        raise RuntimeError("The custom sequences file is empty.")

    if len(custom_fasta_paths.columns) == 1:
        raise RuntimeError(
            "The file you have provided is not TAB delimited. " "Please provide a file with the correct delimiters."
        )

    if 1 < len(custom_fasta_paths.columns) < 3:
        raise RuntimeError(
            "The file you have provided might be missing one of the required fields. "
            "Please provide a file with the correct delimiters, "
            "and the correct number of required fields."
        )

    if len(custom_fasta_paths.columns) > 3:
        raise RuntimeError(
            "The file you have provided might be having more fields than the ones required. "
            "Please provide a file with the correct delimiters, "
            "and the correct number of required fields."
        )

    if custom_fasta_paths["species"].duplicated().any():
        # TODO add a --force flag that allows for this to resolved by haystack
        raise RuntimeError(
            "You have provided more than one sequence for a taxon. "
            "Only one sequence per taxon is allowed. "
            "Please only provide your favourite sequence for each taxon."
        )

    check_unique_taxa_in_custom_input(config["accessions"], config["sequences"])

    inputs = []

    for key, seq in custom_fasta_paths.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["accession"],
        )
        inputs.append(
            config["cache"]
            + "/ncbi/{orgname}/custom_seq-{accession}.fasta.gz".format(orgname=orgname, accession=accession)
        )

    return inputs


rule entrez_aggregate_custom_seqs:
    input:
        get_paths_for_custom_seqs(),
    log:
        config["db_output"] + "/bowtie/custom_seqs.log",
    output:
        config["db_output"] + "/bowtie/custom_seqs.fasta.gz",
    message:
        "Concatenating all the user provided sequences."
    script:
        "../scripts/bowtie2_multifasta.py"


def get_paths_for_custom_acc(wildcards):
    if config["accessions"] == "":
        return ""

    custom_accessions = pd.read_csv(config["accessions"], sep="\t", header=None, names=["species", "accession"])

    if len(custom_accessions) == 0:
        raise RuntimeError("The custom accessions file is empty.")

    if len(custom_accessions.columns) == 1:
        raise RuntimeError(
            "The file you have provided is either not TAB delimited "
            "or it is missing one of the required fields. "
            "Please provide a file with the correct delimiters, "
            "and the correct number of required fields."
        )

    if len(custom_accessions.columns) > 2:
        raise RuntimeError(
            "The file you have provided might be having more fields than the ones required. "
            "Please provide a file with the correct delimiters, "
            "and the correct number of required fields."
        )

    # TODO add a --force flag that allows for this to resolved by haystack
    if custom_accessions["species"].duplicated().any():
        raise RuntimeError(
            "You have provided more than one sequence for a taxon. "
            "Only one sequence per taxon is allowed. "
            "Please only provide your favourite sequence for each taxon."
        )

    check_unique_taxa_in_custom_input(config["accessions"], config["sequences"])

    inputs = []

    for key, seq in custom_accessions.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["accession"],
        )
        inputs.append(
            config["cache"] + "/ncbi/{orgname}/{accession}.fasta.gz".format(orgname=orgname, accession=accession)
        )

    return inputs


rule entrez_aggregate_custom_acc:
    input:
        get_paths_for_custom_acc,
    log:
        config["db_output"] + "/bowtie/custom_acc.log",
    output:
        config["db_output"] + "/bowtie/custom_acc.fasta.gz",
    message:
        "Concatenating all the sequences from user provided accessions."
    script:
        "../scripts/bowtie2_multifasta.py"
