#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import pysam


def mutation_type(alleles):
    """
    Is this mutation a transition (A <-> G and C <-> T) or a transversion (everything else).
    """
    return 'ts' if set(alleles) == {'A', 'G'} or set(alleles) == {'C', 'T'} else 'tv'


def count_ts_tv_init(bam_file, output_file, taxon):
    bam = pysam.AlignmentFile(bam_file, 'rb')

    for read in bam.fetch():
        ts, tv = 0, 0
        for base_call, base_ref in zip(read.seq, read.get_reference_sequence()):
            if base_call != base_ref:
                if mutation_type([base_call, base_ref]) == 'ts':
                    ts += 1
                else:
                    tv += 1

        with open(output_file, 'a') as fout:
            print(taxon, read.query_name, ts, tv, file=fout, sep=",")


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    count_ts_tv_init(
        bam_file=snakemake.input[0],
        output_file=snakemake.output[0],
        taxon=snakemake.wildcards.orgname
    )
