#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd
import sys


def entrez_db_list(acc_input_list, taxa_input_list, output):
    """Function that stores all the taxa and accession is our database in a file"""

    input_list = [[tax, acc] for tax, acc in zip(acc_input_list, taxa_input_list)]
    db_list_df = pd.DataFrame(input_list)
    db_list_df.to_csv(output, sep="\t", header=False, index=False)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    entrez_refseq_create_files(
        acc_input_list=snakemake.input[0], taxa_input_list=snakemake.input[1], output=snakemake.output[0],
    )
