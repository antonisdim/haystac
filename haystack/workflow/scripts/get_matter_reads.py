#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import argparse
import gzip
import sys

import pandas as pd
from Bio import SeqIO


def get_grey_matter_reads(input_fastq, matrix_file, output_fastq, matter):
    """Function to extract all the grey matter reads into a fastq file."""

    dirichlet_matrix = pd.read_csv(matrix_file, sep=",")

    ts_tv_group = dirichlet_matrix.groupby("Read_ID").sum().squeeze()

    if matter == "grey":
        grey_matter_df = ts_tv_group.where(ts_tv_group == 0)
        grey_matter_reads = grey_matter_df[grey_matter_df["Dirichlet_Assignment"].notna()].index.tolist()

        with gzip.open(input_fastq, "rt") as input_handle:
            with gzip.open(output_fastq, "wt") as output_handle:

                for record in SeqIO.parse(input_handle, "fastq"):

                    if record.id in grey_matter_reads:
                        SeqIO.write(record, output_handle, "fastq")

    elif matter == "dark":
        aligned_read_names = ts_tv_group.index.tolist()
        with gzip.open(input_fastq, "rt") as input_handle:
            with gzip.open(output_fastq, "wt") as output_handle:

                for record in SeqIO.parse(input_handle, "fastq"):

                    if record.id not in aligned_read_names:
                        SeqIO.write(record, output_handle, "fastq")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get the grey or dark matter reads..")

    parser.add_argument(
        "--input_fastq", help="Input fastq file to extract reads from.", metavar="",
    )
    parser.add_argument(
        "--matrix_file",
        help="Matrix file produced by haystac to that contains the posterior probabilities of each read belonging to "
        "a taxon in the db.",
        metavar="",
    )

    parser.add_argument(
        "--output_fastq", help="Output file path.", metavar="",
    )

    parser.add_argument(
        "--matter", help="Dark or grey matter reads.", metavar="",
    )

    args = parser.parse_args()

    # noinspection PyUnresolvedReferences
    get_grey_matter_reads(
        input_fastq=args.input_fastq, matrix_file=args.matrix_file, output_fastq=args.output_fastq, matter=args.matter
    )
