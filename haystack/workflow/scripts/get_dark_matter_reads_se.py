#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


from Bio import SeqIO
import pandas as pd
import gzip


def get_dark_matter_reads_se(input_fastq, matrix_file, output_fastq):
    """Function to extract all the grey matter reads into a fastq file."""

    dirichlet_matrix = pd.read_csv(matrix_file, sep=",")

    ts_tv_group = dirichlet_matrix.groupby("Read_ID").sum().squeeze()
    aligned_read_names = ts_tv_group.index.tolist()

    with gzip.open(input_fastq, "rt") as input_handle:
        with gzip.open(output_fastq, "wt") as output_handle:

            for record in SeqIO.parse(input_handle, "fastq"):

                if record.id not in aligned_read_names:
                    SeqIO.write(record, output_handle, "fastq")


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    get_dark_matter_reads_se(
        input_fastq=snakemake.input[0], matrix_file=snakemake.input[1], output_fastq=snakemake.output[0],
    )
