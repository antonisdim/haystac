#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import shutil
import sys


def bowtie_multifasta(fasta_files, output_file):

    print("We're creating the bowtie2 multifasta file for the filtering ...", file=sys.stderr)

    with open(output_file, 'wb') as fout:
        for fasta_file in fasta_files:
            with open(fasta_file, 'rb') as fin:
                shutil.copyfileobj(fin, fout)

    print("Mulitfasta file created", file=sys.stderr)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    bowtie_multifasta(
        fasta_files=snakemake.input,
        output_file=snakemake.output[0]
    )
