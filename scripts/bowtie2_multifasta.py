#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import sys
import gzip
from Bio import bgzf


def bowtie_multifasta(fasta_files, output_file):
    print(
        "We're creating the bowtie2 multifasta file for the filtering ...",
        file=sys.stderr,
    )

    with bgzf.open(output_file, "wt") as fout:
        for fasta_file in fasta_files:
            print(fasta_file, file=sys.stderr)
            with bgzf.open(fasta_file, "rb") as fin:
                shutil.copyfileobj(fin, fout)

    print("Multifasta file created", file=sys.stderr)


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    bowtie_multifasta(fasta_files=snakemake.input, output_file=snakemake.output[0])
