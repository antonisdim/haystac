#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import gzip
import sys

import pandas as pd
from Bio import SeqIO


def get_grey_matter_reads(input_fastq_r1, input_fastq_r2, matrix_file, output_fastq_r1, output_fastq_r2):
    """Function to extract all the grey matter reads into 2 fastq files."""

    dirichlet_matrix = pd.read_csv(matrix_file, sep=",")

    ts_tv_group = dirichlet_matrix.groupby("Read_ID").sum().squeeze()
    grey_matter_df = ts_tv_group.where(ts_tv_group == 0)
    grey_matter_reads = grey_matter_df[grey_matter_df["Dirichlet_Assignment"].notna()].index.tolist()

    with gzip.open(input_fastq_r1, "rt") as input_handle_r1:
        with gzip.open(output_fastq_r1, "wt") as output_handle_r1:

            for record in SeqIO.parse(input_handle_r1, "fastq"):

                if record.id in grey_matter_reads:
                    SeqIO.write(record, output_handle_r1, "fastq")

    with gzip.open(input_fastq_r2, "rt") as input_handle_r2:
        with gzip.open(output_fastq_r2, "wt") as output_handle_r2:

            for record in SeqIO.parse(input_handle_r2, "fastq"):

                if record.id in grey_matter_reads:
                    SeqIO.write(record, output_handle_r2, "fastq")


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    get_grey_matter_reads(
        input_fastq_r1=snakemake.input[0],
        input_fastq_r2=snakemake.input[1],
        matrix_file=snakemake.input[1],
        output_fastq_r1=snakemake.output[0],
        output_fastq_r2=snakemake.output[1],
    )
