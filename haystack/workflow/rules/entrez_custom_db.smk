#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import re

import pandas as pd

from haystack.workflow.scripts.utilities import (
    normalise_name,
    check_unique_taxa_in_custom_input,
    FAIL,
    END,
    is_tty,
)


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
        if not config["resolve_accessions"]:
            dup_taxa = ", ".join(
                [i for i in custom_fasta_paths[custom_fasta_paths["species"].duplicated()]["species"].to_list()]
            )
            err_message = (
                f"You have provided more than one sequence for {dup_taxa}. "
                f"Only one sequence per taxon is allowed. "
                f"Please only provide your favourite sequence for each taxon."
            )
            raise RuntimeError(f"{FAIL}{err_message}{END}" if is_tty else f"{err_message}")
        else:
            custom_fasta_paths = custom_fasta_paths[~custom_fasta_paths["species"].duplicated()]

    check_unique_taxa_in_custom_input(config["accessions"], config["sequences"])

    if config["exclude_accessions"]:
        custom_fasta_paths = custom_fasta_paths[~custom_fasta_paths["accession"].isin(config["exclude_accessions"])]

    inputs = []

    for key, seq in custom_fasta_paths.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            str(seq["accession"]).strip(),
        )

        # TODO make wildcard_constraints consistent with this
        if not re.match("^[\w.]+$", accession):
            err_message = f"The accession '{accession}' for '{orgname}' contains an illegal character"
            raise RuntimeError(f"{FAIL}{err_message}{END}" if is_tty else f"{err_message}")

        inputs.append(config["cache"] + f"/ncbi/{orgname}/custom_seq-{accession}.fasta.gz")

    return inputs


def get_paths_for_custom_acc(wildcards):
    if config["accessions"] == "":
        return ""

    custom_accessions = pd.read_csv(config["accessions"], sep="\t", header=None, names=["species", "accession"])

    if custom_accessions["species"].duplicated().any():
        if not config["resolve_accessions"]:
            dup_taxa = ", ".join(
                [i for i in custom_accessions[custom_accessions["species"].duplicated()]["species"].to_list()]
            )
            err_message = (
                f"You have provided more than one sequence for {dup_taxa}. "
                f"Only one sequence per taxon is allowed. "
                f"Please only provide your favourite sequence for each taxon."
            )
            raise RuntimeError(f"{FAIL}{err_message}{END}" if is_tty else f"{err_message}")
        else:
            custom_accessions = custom_accessions[~custom_accessions["species"].duplicated()]

    check_unique_taxa_in_custom_input(config["accessions"], config["sequences"])

    if config["exclude_accessions"]:
        custom_accessions = custom_accessions[~custom_accessions["accession"].isin(config["exclude_accessions"])]

    inputs = []

    for key, seq in custom_accessions.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            str(seq["accession"]).strip(),
        )

        # TODO make wildcard_constraints consistent with this
        if not re.match("^[\w.]+$", accession):
            err_message = f"The accession '{accession}' for '{orgname}' contains an illegal character"
            raise RuntimeError(f"{FAIL}{err_message}{END}" if is_tty else f"{err_message}")

        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz")

    return inputs
