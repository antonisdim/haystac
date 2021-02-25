#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import beta, hmean


def calculate_dirichlet_abundances(ts_tv_file, p_values_file, total_fastq_reads, sample_abundance):
    """
    Function that calculates the mean posterior abundances of species in metagenomic samples/libraries.
    """

    assert os.stat(ts_tv_file).st_size, f"The ts_tv count file is empty {ts_tv_file}"
    assert os.stat(p_values_file).st_size, f"The chi-square p values file is empty {p_values_file}"
    assert os.stat(total_fastq_reads).st_size, f"The total fastq reads file is empty {total_fastq_reads}"

    # I calculate the coverage, fraction of the covered genome and evenness of coverage of each taxon
    # from reads in its bam/pileup file. Let's go there
    cov_val = pd.read_csv(
        p_values_file,
        sep="\t",
        names=["species", "ref_bases_cov", "total_bases_cov", "coverage", "fraction_ref_cov", "cov_evenness"],
        usecols=["species", "coverage", "fraction_ref_cov", "cov_evenness"],
    )

    evenness_vector = (
        cov_val[["species", "cov_evenness"]]
        .fillna(value=0)
        .groupby("species")
        .apply(hmean)
        .astype("float64")
        .rename("Taxon")
    )
    evenness_vector["Dark_Matter"] = np.nan
    evenness_vector["Grey_Matter"] = np.nan

    fraction_vector = (
        cov_val[["species", "fraction_ref_cov"]]
        .fillna(value=0)
        .groupby("species")
        .apply(hmean)
        .astype("float64")
        .rename("Taxon")
    )
    fraction_vector["Dark_Matter"] = np.nan
    fraction_vector["Grey_Matter"] = np.nan

    coverage_vector = (
        cov_val[["species", "coverage"]]
        .fillna(value=0)
        .groupby("species")
        .apply(hmean)
        .astype("float64")
        .rename("Taxon")
    )
    coverage_vector["Dark_Matter"] = np.nan
    coverage_vector["Grey_Matter"] = np.nan

    ts_tv_matrix = pd.read_csv(ts_tv_file, sep=",", usecols=["Taxon", "Read_ID", "Dirichlet_Assignment"])

    aln_reads_vector = ts_tv_matrix[["Taxon", "Read_ID"]].groupby("Taxon").count().squeeze(axis=1).rename("Taxon")
    aln_reads_vector["Dark_Matter"] = 0
    aln_reads_vector["Grey_Matter"] = 0

    # Sum the Dirichlet Assignments per taxon and calculate the Dark Matter reads from the Dirichlet Assignment column
    ts_tv_group = ts_tv_matrix.groupby("Read_ID").sum().squeeze(axis=1)
    grey_matter = ts_tv_group.where(ts_tv_group == 0).replace(0, 1).fillna(0)

    if len(ts_tv_matrix.Taxon.unique()) > 1:
        a = ts_tv_matrix.groupby("Taxon").sum().squeeze().astype(float)
    else:
        a = ts_tv_matrix.groupby("Taxon").sum().iloc[:, 0].astype(float)
    a.loc["Grey_Matter"] = grey_matter.sum()

    # Add the non aligned filtered reads count in the Dark Matter category
    total_fastq_reads = float(open(total_fastq_reads, "r").read())
    reads_in_bams = len(ts_tv_matrix["Read_ID"].unique())

    remaining_dark_matter = total_fastq_reads - reads_in_bams

    a.loc["Dark_Matter"] = remaining_dark_matter

    print(a, file=sys.stderr)

    # Perform Alberto's formulas
    b = a.sum()

    posterior_abundance_mean = a.add(1).divide(b + len(a)).sort_values(ascending=False)

    # Prepare the dataframe that is going to be outputted and calculate the rest of the output columns.
    posterior_abundance = posterior_abundance_mean.to_frame().reset_index()
    posterior_abundance.rename(columns={"Dirichlet_Assignment": "Mean_Posterior_Abundance"}, inplace=True)
    posterior_abundance["95_CI_lower"] = np.nan
    posterior_abundance["95_CI_upper"] = np.nan
    posterior_abundance["Minimum_Read_Num"] = np.nan
    posterior_abundance["Maximum_Read_Num"] = np.nan
    posterior_abundance["Dirichlet_Read_Num"] = np.nan
    posterior_abundance["Evenness_of_Coverage_Ratio"] = np.nan
    posterior_abundance["Fraction_of_Genome_Covered"] = np.nan
    posterior_abundance["Coverage"] = np.nan
    posterior_abundance["Aligned_Read_Num"] = np.nan

    for idx, row in posterior_abundance.iterrows():
        ai = a.loc[posterior_abundance.iloc[idx, 0]]

        print(ai, file=sys.stderr)

        ci = beta.interval(0.95, ai + 1, b + len(a) - ai - 1)

        print(len(a), file=sys.stderr)

        posterior_abundance.iloc[idx, 2] = ci[0]
        posterior_abundance.iloc[idx, 3] = ci[1]
        posterior_abundance.iloc[idx, 4] = round(ci[0] * b)
        posterior_abundance.iloc[idx, 5] = round(ci[1] * b)
        posterior_abundance.iloc[idx, 6] = a.loc[posterior_abundance.iloc[idx, 0]]
        posterior_abundance.iloc[idx, 7] = evenness_vector.loc[posterior_abundance.iloc[idx, 0]]
        posterior_abundance.iloc[idx, 8] = fraction_vector.loc[posterior_abundance.iloc[idx, 0]]
        posterior_abundance.iloc[idx, 9] = coverage_vector.loc[posterior_abundance.iloc[idx, 0]]
        posterior_abundance.iloc[idx, 10] = aln_reads_vector.loc[posterior_abundance.iloc[idx, 0]]

    with open(sample_abundance, "w") as output_handle:
        posterior_abundance.to_csv(path_or_buf=output_handle, sep="\t", index=False, header=True)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    calculate_dirichlet_abundances(
        ts_tv_file=snakemake.input[0],
        p_values_file=snakemake.input[1],
        total_fastq_reads=snakemake.input[2],
        sample_abundance=snakemake.output[0],
    )
