#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import re

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

    if custom_fasta_paths["species"].duplicated().any():
        # TODO add a --force flag that allows for this to be resolved by haystack
        # TODO tell the user which taxa!
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
            str(seq["accession"]).strip(),
        )

        # TODO make wildcard_constraints consistent with this
        if not re.match("^[\w.]+$", accession):
            raise RuntimeError(f"The accession '{accession}' for '{orgname}' contains an illegal character")

        inputs.append(config["cache"] + f"/ncbi/{orgname}/custom_seq-{accession}.fasta.gz")

    return inputs


rule entrez_sequence_file_aggregator:
    input:
        get_paths_for_custom_seqs(),
    log:
        config["db_output"] + "/bowtie/custom_seqs.log",
    output:
        config["db_output"] + "/bowtie/custom_seqs.done",
    message:
        "Aggregating all the user provided sequences."
    shell:
        "touch {output}"


def get_paths_for_custom_acc(wildcards):
    if config["accessions"] == "":
        return ""

    custom_accessions = pd.read_csv(config["accessions"], sep="\t", header=None, names=["species", "accession"])

    # TODO add a --force flag that allows for this to be resolved by haystack
    # TODO tell the user which taxa!
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
            str(seq["accession"]).strip(),
        )

        # TODO make wildcard_constraints consistent with this
        if not re.match("^[\w.]+$", accession):
            raise RuntimeError(f"The accession '{accession}' for '{orgname}' contains an illegal character")

        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz")

    return inputs


rule entrez_accessions_file_aggregator:
    input:
        get_paths_for_custom_acc,
    log:
        config["db_output"] + "/bowtie/custom_acc.log",
    output:
        config["db_output"] + "/bowtie/custom_acc.done",
    message:
        "Aggregating all the sequences from user provided accessions."
    shell:
        "touch {output}"
