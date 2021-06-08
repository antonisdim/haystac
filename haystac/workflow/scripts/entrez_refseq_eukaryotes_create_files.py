#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import sys

import numpy as np
import pandas as pd

from haystac.workflow.scripts.utilities import REGEX_BLACKLIST


def entrez_refseq_eukaryotes_create_files(
    config,
    input_file,
    euk_genomes_out,
):

    """Function to parse the refseq genomes report for eukaryotes."""

    # read the file

    refseq_eukaryotes = pd.read_csv(input_file, sep="\t")

    # drop duplicate species/strains and get pairs of taxa and a assembly accession codes

    refseq_eukaryotes_rmdup = refseq_eukaryotes[~refseq_eukaryotes["#Organism/Name"].duplicated(keep="first")]

    eukaryotes_unique = refseq_eukaryotes_rmdup[["#Organism/Name", "Assembly Accession"]].copy()

    # drop rows that have no accessions

    means_empty_record = ["", "-", np.nan]

    eukaryotes = eukaryotes_unique.loc[~eukaryotes_unique["Assembly Accession"].isin(means_empty_record)]

    # rename columns

    eukaryotes.rename(columns={"#Organism/Name": "species", "Assembly Accession": "AccessionVersion"}, inplace=True)

    # regex for species name

    eukaryotes["species"] = eukaryotes["species"].replace(REGEX_BLACKLIST, "_", regex=True)

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
            eukaryotes = eukaryotes[(~eukaryotes["species"].isin(user_df["species"]))]

    # print the output to csv

    header = ["species", "AccessionVersion"]
    eukaryotes.to_csv(euk_genomes_out, sep="\t", header=header, index=False)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    entrez_refseq_eukaryotes_create_files(
        config=snakemake.config,
        input_file=snakemake.input[0],
        euk_genomes_out=snakemake.output[0],
    )
