#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import sys


def calculate_bt2_idx_chunks(mem_resources, mem_rescaling_factor, fasta_files, output):
    """Calculate the number of chunks that the db sequences are going to be split into"""

    chunk_size = float(mem_resources) / float(mem_rescaling_factor)
    total_size = 0.0
    chunks = 1

    # TODO handle edge case where mem_resources is so low we can't build an index for a single large fasta file

    for fasta_file in fasta_files:
        file_size = os.stat(fasta_file).st_size / (1024 ** 2)

        if total_size + file_size >= chunk_size:
            total_size = file_size
            chunks += 1
        else:
            total_size += file_size

    with open(output, "w") as outfile:
        print(chunks, file=outfile)


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    calculate_bt2_idx_chunks(
        mem_resources=snakemake.params[1],
        mem_rescaling_factor=snakemake.params[2],
        fasta_files=snakemake.input,
        output=snakemake.output[0],
    )
