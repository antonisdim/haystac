#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import json
import os
import sys

import numpy as np
import pandas as pd

from math import ceil
from haystac.workflow.scripts.utilities import PE


def calculate_likelihoods(ts_tv_file, readlen_file, taxa_file_paths, config, output_matrix, output_params):
    """
    Calculate the parameters for the analytical framework of the method, then the likelihoods for each read/taxon pair
    and perform the dirichlet distribution assignment.
    """
    assert os.stat(ts_tv_file).st_size, f"The ts_tv count file is empty {ts_tv_file}"
    assert os.stat(readlen_file).st_size, f"The read length is empty {readlen_file}"
    assert len(taxa_file_paths) > 0, f"The taxa list is empty {taxa_file_paths}"

    print("Reading the initial Ts/Tv matrix.", file=sys.stderr)
    init_ts_tv = pd.read_csv(ts_tv_file, names=["Taxon", "Read_ID", "Ts", "Tv"], sep=",")

    print("Calculating the sum of the transitions.", file=sys.stderr)
    ts_sum = int(init_ts_tv["Ts"].sum())

    print("Calculating the sum of the transversions.", file=sys.stderr)
    tv_sum = int(init_ts_tv["Tv"].sum())

    aligned_read_count = len(init_ts_tv["Read_ID"].unique())
    average_read_length = float(open(readlen_file, "r").read())

    max_mismatch = ceil(config["mismatch_probability"] * float(average_read_length))

    # calculate the parameters for the analytical framework
    if ts_sum == 0 and tv_sum != 0:
        sigma_t = 0
        sigma_v = config["mismatch_probability"]

        delta_t = 0
        delta_v = sigma_v / float(1 - sigma_v - sigma_t)

        ts_missing_val = 0
        tv_missing_val = max_mismatch

    elif ts_sum != 0 and tv_sum == 0:
        sigma_t = config["mismatch_probability"]
        sigma_v = 0

        delta_t = sigma_t / float(1 - sigma_v - sigma_t)
        delta_v = 0

        ts_missing_val = max_mismatch
        tv_missing_val = 0

    elif ts_sum == 0 and tv_sum == 0:
        sigma_t = config["mismatch_probability"] / (float(2))
        sigma_v = config["mismatch_probability"] / (float(2))

        delta_v = sigma_v / float(1 - sigma_v - sigma_t)
        delta_t = sigma_t / float(1 - sigma_v - sigma_t)

        ts_missing_val = round(max_mismatch / float(2))
        tv_missing_val = round(max_mismatch / float(2))

    else:
        ts_tv_ratio = ts_sum / float(tv_sum)
        fixed_mismatch_probability = config["mismatch_probability"]

        # solve linear system to find σv and σt
        coef = np.array([[1, 1], [1, -ts_tv_ratio]])
        const = np.array([fixed_mismatch_probability, 0])
        x = np.linalg.solve(coef, const)

        sigma_t = x[0]
        sigma_v = x[1]

        delta_v = sigma_v / float(1 - sigma_v - sigma_t)
        delta_t = sigma_t / float(1 - sigma_v - sigma_t)

        # solve linear system to find missing values for Ts and Tv
        coef_miss_tv = np.array([[1, 1], [1, -ts_tv_ratio]])
        const_miss_tv = np.array([max_mismatch, 0])
        y = np.linalg.solve(coef_miss_tv, const_miss_tv)

        ts_missing_val = round(y[0])
        tv_missing_val = round(y[1])

    # now calculate the likelihoods for each read/taxon pair
    data_ts_missing = pow(delta_t, ts_missing_val)
    data_tv_missing = pow(delta_v, tv_missing_val)

    # as we count max mismatches for every mate of a pair. So if a whole pair is missing we double the missing values
    if config["read_mode"] == PE:
        data_ts_missing = 2 * data_ts_missing
        data_tv_missing = 2 * data_tv_missing

    print("calculating the likelihood numerator", file=sys.stderr)
    init_ts_tv["ll_nom"] = init_ts_tv["Ts"].rpow(delta_t) * init_ts_tv["Tv"].rpow(delta_v)

    total_taxa_count = len(taxa_file_paths)

    print("calculating the proper likelihood", file=sys.stderr)
    init_ts_tv["Likelihood"] = np.nan

    init_ts_tv["Likelihood"] = init_ts_tv.groupby("Read_ID")["ll_nom"].transform(
        lambda nom: nom / (sum(nom) + ((total_taxa_count - nom.count()) * data_ts_missing * data_tv_missing))
    )

    print(init_ts_tv, file=sys.stderr)

    # every likelihood equal or above 0.95 gets assigned as 1, everything below as 0
    init_ts_tv["Dirichlet_Assignment"] = (
        init_ts_tv["Likelihood"]
        .where(init_ts_tv["Likelihood"] > config["min_prob"], 0.0)
        .where(init_ts_tv["Likelihood"] <= config["min_prob"], 1.0)
    )

    init_ts_tv.to_csv(output_matrix, sep=",", index=False)

    with open(output_params, "w") as fout:
        params = {
            "T": ts_sum,
            "V": tv_sum,
            "sigma_v": sigma_v,
            "sigma_t": sigma_t,
            "delta_v": delta_v,
            "delta_t": delta_t,
            "ts_missing_val": ts_missing_val,
            "tv_missing_val": tv_missing_val,
            "aligned_read_count": aligned_read_count,
        }

        json.dump(params, fout)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    sys.stderr = open(snakemake.log[0], "w")

    # noinspection PyUnresolvedReferences
    calculate_likelihoods(
        ts_tv_file=snakemake.input[0],
        readlen_file=snakemake.input[1],
        taxa_file_paths=snakemake.input[2],
        config=snakemake.config,
        output_matrix=snakemake.output[0],
        output_params=snakemake.output[1],
    )
