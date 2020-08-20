#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import sys
import gzip
import pandas as pd
import shutil

from Bio import bgzf


# script finds the input file

# reads the file


def entrez_custom_sequences(config, taxon, output_file):
    if os.path.exists(config["sequences"]):

        custom_fasta_paths = pd.read_csv(
            config["sequences"], sep="\t", header=None, names=["species", "accession", "path"],
        )

    else:
        raise RuntimeError(
            "The file containing the paths to the custom fasta sequences isn not there. "
            "Please provide a valid path for the tab delimited input file."
        )

    fasta_file = custom_fasta_paths.loc[custom_fasta_paths["species"] == taxon]["path"].values[0]

    if os.path.exists(fasta_file):
        print("File for taxon {taxon} exists.".format(taxon=taxon), file=sys.stderr)
        filename, file_extension = os.path.splitext(fasta_file)
        if file_extension == ".gz":
            print(
                "We're putting the fasta file for taxon {taxon} in the database.".format(taxon=taxon), file=sys.stderr,
            )

            with bgzf.open(output_file, "wt") as fout:
                print(fasta_file, file=sys.stderr)
                with gzip.open(fasta_file, "rb") as fin:
                    shutil.copyfileobj(fin, fout)

            print("Compressed fasta file created", file=sys.stderr)

        else:
            print(
                "We're putting the fasta file for taxon {taxon} in the database.".format(taxon=taxon), file=sys.stderr,
            )

            with bgzf.open(output_file, "wt") as fout:
                print(fasta_file, file=sys.stderr)
                with open(fasta_file, "r") as fin:
                    shutil.copyfileobj(fin, fout)

    else:
        raise RuntimeError(
            "The path for the fasta file input for taxon {taxon} isn't valid. "
            "PLease provide a valid path in your input file.".format(taxon=taxon)
        )


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    entrez_custom_sequences(
        config=snakemake.config, taxon=snakemake.wildcards.orgname, output_file=snakemake.output[0],
    )
