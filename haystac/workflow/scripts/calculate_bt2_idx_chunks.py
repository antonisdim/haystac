#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os

import pandas as pd

from haystac.workflow.scripts.utilities import print_error


def calculate_bt2_idx_chunks(mem_resources, mem_rescale_factor, fasta_files, output_tsv, output_txt):
    """Calculate the number of chunks that the db sequences are going to be split into"""

    fasta_paths_random = []
    with open(str(fasta_files), "r") as fin:
        for line in fin:
            fasta_paths_random.append(line.strip())

    chunk_size = float(mem_resources) / float(mem_rescale_factor)
    total_size = 0.0
    chunks = 1

    chunk_df_list = []

    for fasta_file in fasta_paths_random:
        file_size = os.stat(fasta_file).st_size / (1024 ** 2)

        if file_size >= mem_resources or file_size >= chunk_size:
            print_error(
                f"Fasta file {fasta_file} is bigger than the RAM resources provided. "
                f"Unfortunately an index cannot be built."
            )

        if total_size + file_size >= chunk_size:
            total_size = file_size
            chunks += 1
        else:
            total_size += file_size

        chunk_df_list.append([chunks, fasta_file])

    chunk_df = pd.DataFrame(chunk_df_list, columns=["chunk", "path"])

    chunk_df.to_csv(output_tsv, sep="\t", index=False, header=False)

    idx_chunk_total = chunk_df["chunk"].max()

    with open(output_txt, "w") as outfile:
        print(idx_chunk_total, file=outfile)


if __name__ == "__main__":

    # noinspection PyUnresolvedReferences
    calculate_bt2_idx_chunks(
        mem_resources=snakemake.params[1],
        mem_rescale_factor=snakemake.params[2],
        fasta_files=snakemake.input,
        output_tsv=snakemake.output[0],
        output_txt=snakemake.output[1],
    )
