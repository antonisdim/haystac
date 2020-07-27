#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

##### Target rules #####


from scripts.rip_utilities import normalise_name, check_unique_taxa_in_custom_input


rule entrez_custom_sequences:
    input:
        config["sequences"],
    log:
        config["genome_cache_folder"] + "/{orgname}/custom_seq-{accession}.log",
    output:
        config["genome_cache_folder"] + "/{orgname}/custom_seq-{accession}.fasta.gz",
    message:
        "Incorporating the user provided fasta sequence {wildcards.accession} for taxon {wildcards.orgname}. "
        "The provided sequence can be found in {output} and its log file in {log}."
    script:
        "../scripts/entrez_custom_sequences.py"


def get_paths_for_custom_seqs():

    custom_fasta_paths = pd.read_csv(
        config["sequences"],
        sep="\t",
        header=None,
        names=["species", "accession", "path"],
    )

    if len(custom_fasta_paths) == 0:
        raise RuntimeError("The custom sequences file is empty.")

    if custom_fasta_paths["species"].duplicated().any():
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
            config["genome_cache_folder"]
            + "/{orgname}/custom_seq-{accession}.fasta.gz".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


rule entrez_aggregate_custom_seqs:
    input:
        get_paths_for_custom_seqs,
    log:
        config["db_output"] + "/bowtie/{query}_custom_seqs.log",
    output:
        config["db_output"] + "/bowtie/{query}_custom_seqs.fasta.gz",
    message:
        "Concatenating all the user provided sequences {input} in {output}. "
        "Its log file can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"


def get_paths_for_custom_acc():

    custom_accessions = pd.read_csv(
        config["custom_acc_file"], sep="\t", header=None, names=["species", "accession"]
    )

    if len(custom_accessions) == 0:
        raise RuntimeError("The custom accessions file is empty.")

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
            config["genome_cache_folder"]
            + "/{orgname}/{accession}.fasta.gz".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


rule entrez_aggregate_custom_acc:
    input:
        get_paths_for_custom_acc,
    log:
        config["db_output"] + "/bowtie/{query}_custom_acc.log",
    output:
        config["db_output"] + "/bowtie/{query}_custom_acc.fasta.gz",
    message:
        "Concatenating all the sequences from user provided accessions {input} in {output}. "
        "Its log file can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"
