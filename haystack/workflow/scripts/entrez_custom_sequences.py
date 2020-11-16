#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import gzip
import os
import shutil
import sys

import pandas as pd
from Bio import bgzf

from haystack.workflow.scripts.utilities import print_error


def entrez_custom_sequences(config, taxon, output_file):
    if os.path.exists(config["sequences"]):

        custom_fasta_paths = pd.read_csv(
            config["sequences"], sep="\t", header=None, names=["species", "accession", "path"],
        )

    else:
        print_error(
            "The file containing the paths to the custom fasta sequences isn not there. "
            "Please provide a valid path for the tab delimited input file."
        )

    # noinspection PyUnboundLocalVariable
    fasta_file = custom_fasta_paths.loc[custom_fasta_paths["species"] == taxon]["path"].values[0]

    if os.path.exists(fasta_file):
        print(f"File for taxon {taxon} exists.", file=sys.stderr)
        filename, file_extension = os.path.splitext(fasta_file)
        if file_extension == ".gz":
            print(
                f"We're putting the fasta file for taxon {taxon} in the database.", file=sys.stderr,
            )

            with bgzf.open(output_file, "wt") as fout:
                print(fasta_file, file=sys.stderr)
                with gzip.open(fasta_file, "rb") as fin:
                    shutil.copyfileobj(fin, fout)

            print("Compressed fasta file created", file=sys.stderr)

        else:
            print(
                f"We're putting the fasta file for taxon {taxon} in the database.", file=sys.stderr,
            )

            with bgzf.open(output_file, "wt") as fout:
                print(fasta_file, file=sys.stderr)
                with open(fasta_file, "r") as fin:
                    shutil.copyfileobj(fin, fout)

    else:
        print_error(
            f"The path for the fasta file input for taxon {taxon} "
            f"isn't valid. Please provide a valid path in your input file."
        )


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    entrez_custom_sequences(
        config=snakemake.config, taxon=snakemake.wildcards.orgname, output_file=snakemake.output[0],
    )
