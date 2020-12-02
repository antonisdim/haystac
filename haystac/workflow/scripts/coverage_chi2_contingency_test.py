#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import sys

import pandas as pd
from scipy.stats import chi2_contingency


def coverage_chi2_contingency_test(cov_file, taxon_fasta_idx, taxon, outfile):
    """
    Function that calculates the pvalue of the coverage. Testing if there was clustering bias during the
    sequencing/the reads that contribute to the identification/abundance of a species come only
    from a specific genomic region
    """

    assert os.stat(cov_file).st_size, f"The file with the coverage stats {cov_file} is empty."

    taxon_seqlen = genome_sizes(taxon_fasta_idx)

    cov_stats = pd.read_csv(cov_file, sep="\t", names=["observed", "expected"])

    expected_coverage = cov_stats['expected'].iloc[0]
    observed_coverage = cov_stats['observed'].iloc[0]

    contingency_first_row = [observed_coverage, expected_coverage]
    print("Observed and expected coverage are ", contingency_first_row, file=sys.stderr)

    # oddsratio, pvalue = fisher_exact([contingency_first_row, contingency_second_row])

    chi2, pvalue, dof, expected = chi2_contingency(contingency_first_row)

    with open(outfile, "w") as outhandle:
        print(taxon, pvalue, file=outhandle, sep="\t")


def genome_sizes(taxon_fasta_idx):
    faidx = pd.read_csv(
        taxon_fasta_idx,
        sep="\t",
        names=["Name", "Length", "Offset", "Linebases", "Linewidth"],
    )

    taxon_seq_len = faidx["Length"].sum()

    return taxon_seq_len


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    coverage_chi2_contingency_test(
        cov_file=snakemake.input[0],
        taxon_fasta_idx=snakemake.input[1],
        taxon=snakemake.wildcards.orgname,
        outfile=snakemake.output[0],
    )
