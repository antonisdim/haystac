#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import sys

import pandas as pd
import pysam


def get_dirichlet_reads(input_bam, output_bam, matrix_file, taxon):
    """Function to get dirichlet assigned reads from a taxon's bam file"""

    dirichlet_matrix = pd.read_csv(matrix_file, sep=",")

    taxon_dirichlet_df = dirichlet_matrix[
        ((dirichlet_matrix["Taxon"] == taxon) & (dirichlet_matrix["Dirichlet_Assignment"] == 1))
    ]

    dirichlet_reads = taxon_dirichlet_df["Read_ID"].tolist()

    bam = pysam.AlignmentFile(input_bam, "rb")

    names = set(dirichlet_reads)

    with pysam.AlignmentFile(output_bam, "wb", header=bam.header) as fout:

        for read in bam.fetch():

            if read.query_name in names:
                fout.write(read)


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    get_dirichlet_reads(
        input_bam=snakemake.input[0],
        output_bam=snakemake.output[0],
        matrix_file=snakemake.input[1],
        taxon=snakemake.wildcards.orgname,
    )
