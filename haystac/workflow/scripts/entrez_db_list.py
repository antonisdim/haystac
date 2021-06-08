#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import sys

import pandas as pd


def entrez_db_list(input_list, output):
    """Function that stores all the taxa and accession is our database in a file"""

    db_list = []

    for fasta_file in input_list:
        path_list = fasta_file.split("/")
        taxon = path_list[-2]
        accession = path_list[-1].replace(".fasta.gz", "")
        db_list.append([taxon, accession])

    db_list_df = pd.DataFrame(db_list)
    db_list_df.to_csv(output, sep="\t", header=False, index=False)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    entrez_db_list(
        input_list=snakemake.input,
        output=snakemake.output[0],
    )
