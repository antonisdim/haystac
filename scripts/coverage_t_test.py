#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import pandas as pd
import pysam
from scipy.stats import fisher_exact


def coverage_t_test(bam, selected_seqs_file, nuccore_file, taxon, outfile):
    # TODO add a block comment explaining what this function does

    taxon_seqlen_dict = genome_sizes(selected_seqs_file, nuccore_file)

    # 86) Define function that calculates the pvalue for the coverage
    #     (testing if there was clustering bias during the sequencing)

    pileup_dict = {}

    bamfile = pysam.AlignmentFile(bam, "rb")

    for pileupcolumn in bamfile.pileup():
        pileup_dict[pileupcolumn.pos] = pileupcolumn.n

    expected_coverage = sum(list(pileup_dict.values()))
    observed_coverage = len(list(pileup_dict.keys()))

    contingency_first_row = [observed_coverage, expected_coverage]
    print("Observed and expected coverage are ", contingency_first_row, file=sys.stderr)

    contingency_second_row = [taxon_seqlen_dict[taxon], taxon_seqlen_dict[taxon]]

    oddsratio, pvalue = fisher_exact([contingency_first_row, contingency_second_row])

    with open(outfile, 'w') as outhandle:
        print(taxon, pvalue, file=outhandle, sep='\t')


# TODO this is a really inefficient way of getting the size of the sequences!!
#      we shouldn't have to load either of these files!
def genome_sizes(selected_seqs_file, nuccore_file):
    selected_seqs = pd.read_csv(selected_seqs_file, sep='\t')

    nuccore_seqs = pd.read_csv(nuccore_file, sep='\t')

    merge = pd.merge(selected_seqs, nuccore_seqs, on=['GBSeq_accession-version'])[['species', 'GBSeq_length']]. \
        replace(' ', '_', regex=True)

    taxon_seqlen_dict = pd.Series(merge.GBSeq_length.values, index=merge.species).to_dict()

    return taxon_seqlen_dict


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    coverage_t_test(
        bam=snakemake.input[0],
        selected_seqs_file=snakemake.input[1],
        nuccore_file=snakemake.input[2],
        taxon=snakemake.wildcards.orgname,
        outfile=snakemake.output[0]
    )
