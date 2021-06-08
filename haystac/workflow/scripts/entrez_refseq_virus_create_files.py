#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import sys

import pandas as pd

from haystac.workflow.scripts.utilities import REGEX_BLACKLIST, print_error


def entrez_refseq_virus_create_files(
    config,
    input_file,
    viral_genomes_out,
):

    """Function to parse the refseq genomes report for viruses."""

    # read the file

    refseq_viruses = pd.read_csv(input_file, sep="\t")

    # drop duplicate species/strains

    refseq_viruses_rmdup = refseq_viruses[~refseq_viruses["#Organism/Name"].duplicated(keep="first")]

    viruses_unique = refseq_viruses_rmdup[["#Organism/Name", "Segmemts"]].copy()
    viruses_unique["Accession"] = ""

    # assign a segment accession code to a species.
    # 1 segment acc should be enough as all of them will be fetched through the assembly database.

    for index, row in viruses_unique.iterrows():

        seq_list = row["Segmemts"].split("; ")

        if len(seq_list) == 1 and ":" not in seq_list[0]:
            if "/" in seq_list[0]:
                row["Accession"] = seq_list[0].split("/")[0]
            elif "-" == seq_list[0]:
                continue
            else:
                row["Accession"] = seq_list[0]
        elif len(seq_list) == 1 and ":" in seq_list[0]:
            row["Accession"] = seq_list[0].split(":")[1].split("/")[0]
        elif len(seq_list) > 1:
            if "/" in seq_list[0].split(":")[1]:
                row["Accession"] = seq_list[0].split(":")[1].split("/")[0]
            else:
                row["Accession"] = seq_list[0].split(":")[1]
        else:
            print_error(seq_list)
            break

    # rename columns

    viruses = viruses_unique[["#Organism/Name", "Accession"]]
    viruses.rename(columns={"#Organism/Name": "species", "Accession": "AccessionVersion"}, inplace=True)

    # drop rows that have no accessions

    viruses = viruses[viruses["AccessionVersion"] != ""]

    # regex for species name

    viruses["species"] = viruses["species"].replace(REGEX_BLACKLIST, "_", regex=True)

    # check for duplicates from user input

    if config["sequences"] or config["accessions"]:
        user_inputs = []
        if os.path.isfile(config["sequences"]):
            custom_fasta_paths = pd.read_csv(
                config["sequences"],
                sep="\t",
                header=None,
                names=["species", "accession", "path"],
            )
            user_inputs.append(custom_fasta_paths)
        if os.path.isfile(config["accessions"]):
            custom_accessions = pd.read_csv(
                config["accessions"],
                sep="\t",
                header=None,
                names=["species", "accession"],
            )
            user_inputs.append(custom_accessions)

        for user_df in user_inputs:
            viruses = viruses[(~viruses["species"].isin(user_df["species"]))]

    # print the output to csv

    header = ["species", "AccessionVersion"]
    viruses.to_csv(viral_genomes_out, sep="\t", header=header, index=False)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    entrez_refseq_virus_create_files(
        config=snakemake.config,
        input_file=snakemake.input[0],
        viral_genomes_out=snakemake.output[0],
    )
