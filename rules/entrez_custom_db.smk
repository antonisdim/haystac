#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

WITH_CUSTOM_SEQUENCES = config["WITH_CUSTOM_SEQUENCES"]
WITH_CUSTOM_ACCESSIONS = config["WITH_CUSTOM_ACCESSIONS"]


rule entrez_custom_sequences:
    input:
        config["custom_seq_file"],
    log:
        "database/{orgname}/custom_seq-{accession}.log",
    output:
        "database/{orgname}/custom_seq-{accession}.fasta.gz",
    message:
        "Incorporating the user provided fasta sequence {wildcards.accession} for taxon {wildcards.orgname}. "
        "The provided sequence can be found in {output} and its log file in {log}."
    script:
        "../scripts/entrez_custom_sequences.py"


def get_paths_for_custom_seqs(wildcards):  # TODO wildcards not used

    custom_fasta_paths = pd.read_csv(
        config["custom_seq_file"],
        sep="\t",
        header=None,
        names=["species", "accession", "path"],
    )

    if len(custom_fasta_paths) == 0:
        raise RuntimeError("The custom sequences file is empty.")

    # TODO code duplicated on lines 107-138!
    if custom_fasta_paths["species"].duplicated().any():
        raise RuntimeError(
            "You have provided more than one sequence for a taxon. "
            "Only one sequence per taxon is allowed. "
            "Please only provide your favourite sequence for each taxon."
        )

    if WITH_CUSTOM_ACCESSIONS and WITH_CUSTOM_SEQUENCES:
        custom_fasta_paths = pd.read_csv(
            config["custom_seq_file"],
            sep="\t",
            header=None,
            names=["species", "accession", "path"],
        )
        custom_accessions = pd.read_csv(
            config["custom_acc_file"],
            sep="\t",
            header=None,
            names=["species", "accession"],
        )

        taxon_acc = custom_accessions["species"].tolist()
        taxon_seq = custom_fasta_paths["species"].tolist()

        if bool(set(taxon_acc) & set(taxon_seq)):
            raise RuntimeError(
                "You have provided the same taxon both in your custom sequences file and your "
                "custom accessions file. Please pick and keep ONLY one entry from both of these files. "
                "You can only have 1 sequence per chosen taxon in your database."
            )

    inputs = []

    for key, seq in custom_fasta_paths.iterrows():
        # TODO move the .replace() code into a function called normalize_name() and update any other instances
        orgname, accession = (
            seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
            seq["accession"],
        )
        inputs.append(
            "database/{orgname}/custom_seq-{accession}.fasta.gz".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


rule entrez_aggregate_custom_seqs:
    input:
        get_paths_for_custom_seqs,
    log:
        "{query}/bowtie/{query}_custom_seqs.log",
    output:
        "{query}/bowtie/{query}_custom_seqs.fasta.gz",
    message:
        "Concatenating all the user provided sequences {input} in {output}. "
        "Its log file can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"


def get_paths_for_custom_acc(wildcards):

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

    if WITH_CUSTOM_ACCESSIONS and WITH_CUSTOM_SEQUENCES:
        custom_fasta_paths = pd.read_csv(
            config["custom_seq_file"],
            sep="\t",
            header=None,
            names=["species", "accession", "path"],
        )
        custom_accessions = pd.read_csv(
            config["custom_acc_file"],
            sep="\t",
            header=None,
            names=["species", "accession"],
        )

        taxon_acc = custom_accessions["species"].tolist()
        taxon_seq = custom_fasta_paths["species"].tolist()

        if bool(set(taxon_acc) & set(taxon_seq)):
            raise RuntimeError(
                "You have provided the same taxon both in your custom sequences file and your "
                "custom accessions file. Please pick and keep ONLY one entry from both of these files. "
                "You can only have 1 sequence per chosen taxon in your database."
            )

    inputs = []

    for key, seq in custom_accessions.iterrows():
        # TODO move the .replace() code into a function called normalize_name() and update any other instances
        orgname, accession = (
            seq["species"].replace(" ", "_").replace("[", "").replace("]", ""),
            seq["accession"],
        )
        inputs.append(
            "database/{orgname}/{accession}.fasta.gz".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


rule entrez_aggregate_custom_acc:
    input:
        get_paths_for_custom_acc,
    log:
        "{query}/bowtie/{query}_custom_acc.log",
    output:
        "{query}/bowtie/{query}_custom_acc.fasta.gz",
    message:
        "Concatenating all the sequences from user provided accessions {input} in {output}. "
        "Its log file can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"
