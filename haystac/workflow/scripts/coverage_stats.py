#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os

import pandas as pd


def coverage_chi2_contingency_test(cov_file, taxon_fasta_idx, taxon, outfile):
    """
    Function that calculates the pvalue of the coverage. Testing if there was clustering bias during the
    sequencing/the reads that contribute to the identification/abundance of a species come only
    from a specific genomic region
    """

    assert os.stat(cov_file).st_size, f"The file with the coverage stats {cov_file} is empty."

    cov_stats = pd.read_csv(cov_file, sep="\t", names=["all_positions", "observed", "expected"])

    all_positions = cov_stats["all_positions"].iloc[0]
    ref_bases_cov = cov_stats["observed"].iloc[0]
    total_bases_cov = cov_stats["expected"].iloc[0]

    # calculate coverage
    coverage = total_bases_cov / all_positions

    # calculate the fraction of the genome that is covered
    taxon_seqlen = genome_sizes(taxon_fasta_idx)
    fraction_ref_cov = ref_bases_cov / taxon_seqlen

    # calculate an evenness of coverage metric valid only for data < 1x. If less than 10 / 1 that is good.
    cov_evenness = coverage / fraction_ref_cov

    with open(outfile, "w") as outhandle:
        print(taxon, ref_bases_cov, total_bases_cov, coverage, fraction_ref_cov, cov_evenness, file=outhandle, sep="\t")


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
    coverage_chi2_contingency_test(
        cov_file=snakemake.input[0],
        taxon_fasta_idx=snakemake.input[1],
        taxon=snakemake.wildcards.orgname,
        outfile=snakemake.output[0],
    )
