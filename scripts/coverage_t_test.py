#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import pandas as pd
import pysam
from scipy.stats import fisher_exact


def coverage_t_test(bam, taxon_fasta_idx, taxon, outfile):
    """
    Function that calculates the pvalue of the coverage. Testing if there was clustering bias during the
    sequencing/the reads that contribute to the identification/abundance of a species come only
    from a specific genomic region
    """
    taxon_seqlen = genome_sizes(taxon_fasta_idx)

    pileup_dict = {}

    bamfile = pysam.AlignmentFile(bam, "rb")

    for pileupcolumn in bamfile.pileup():
        pileup_dict[pileupcolumn.pos] = pileupcolumn.n

    expected_coverage = sum(list(pileup_dict.values()))
    observed_coverage = len(list(pileup_dict.keys()))

    contingency_first_row = [observed_coverage, expected_coverage]
    print("Observed and expected coverage are ", contingency_first_row, file=sys.stderr)

    contingency_second_row = [taxon_seqlen, taxon_seqlen]

    oddsratio, pvalue = fisher_exact([contingency_first_row, contingency_second_row])

    with open(outfile, 'w') as outhandle:
        print(taxon, pvalue, file=outhandle, sep='\t')


def genome_sizes(taxon_fasta_idx):
    faidx = pd.read_csv(taxon_fasta_idx, sep='\t', names=['Name', 'Length', 'Offset', 'Linebases', 'Linewidth'])

    taxon_seq_len = faidx['Length'].sum()

    return taxon_seq_len


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    coverage_t_test(
        bam=snakemake.input[0],
        taxon_fasta_idx=snakemake.input[1],
        taxon=snakemake.wildcards.orgname,
        outfile=snakemake.output[0]
    )
