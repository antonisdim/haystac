#!/usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import shutil
import sys


def create_multifasta(outfilename, refseq_folder):

    with open(outfilename, 'wb') as outfile:
        for filename in glob.iglob(refseq_folder + '**/*.fasta', recursive=True):
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    print("We're creating the bowtie2 multifasta file for the filtering ...", file=sys.stderr)

    create_multifasta(
        outfilename=snakemake.output[0],
        refseq_folder=snakemake.input[0]
    )

    print("Mulitfasta file created", file=sys.stderr)
