#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import numpy as np
import pandas as pd
import json
import os


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)


def calculate_likelihoods(ts_tv_file, readlen_file, taxa_file, config, output_matrix, output_params):
    """
    Function whose first part calculates the parameters for the analytical framework of the method.
    The second part calculates the likelihoods for each read/taxon pair and performs the dirichlet distribution
    assignment.

    """

    # First part

    if os.stat(ts_tv_file).st_size == 0:
        raise RuntimeError("The ts_tv count file is empty. Go back and check why.")

    if os.stat(readlen_file).st_size == 0:
        raise RuntimeError("The read length is empty. Go back and check why.")

    if os.stat(taxa_file).st_size == 0:
        raise RuntimeError("The taxa list is empty. Go back and check why.")

    print("Reading the initial Ts/Tv matrix.", file=sys.stderr)
    init_ts_tv = pd.read_csv(ts_tv_file, names=['Taxon', 'Read_ID', 'Ts', 'Tv'], sep=',')

    print('Calculating the sum of the transitions.', file=sys.stderr)
    t = init_ts_tv['Ts'].sum()
    print('Calculating the sum of the transversions.', file=sys.stderr)
    v = init_ts_tv['Tv'].sum()

    aligned_read_count = len(init_ts_tv['Read_ID'].unique())
    average_read_length = float(open(readlen_file, 'r').read())

    max_mismatch = round(config['mismatch_probability'] * float(average_read_length))

    if t == 0 and v != 0:
        sigma_t = 0
        sigma_v = config['mismatch_probability']
        delta_t = 0
        delta_v = sigma_v / float(1 - sigma_v - sigma_t)
        ts_missing_val = 0
        tv_missing_val = max_mismatch
    elif t != 0 and v == 0:
        sigma_t = config['mismatch_probability']
        sigma_v = 0
        delta_t = sigma_t / float(1 - sigma_v - sigma_t)
        delta_v = 0
        ts_missing_val = max_mismatch
        tv_missing_val = 0
    elif t == 0 and v == 0:
        sigma_t = config['mismatch_probability'] / (float(2))
        sigma_v = config['mismatch_probability'] / (float(2))

        delta_v = sigma_v / float(1 - sigma_v - sigma_t)
        delta_t = sigma_t / float(1 - sigma_v - sigma_t)

        ts_missing_val = round(max_mismatch/float(2))
        tv_missing_val = round(max_mismatch/float(2))
    else:
        ts_tv_ratio = t / float(v)
        fixed_mismatch_probability = config['mismatch_probability']
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

    param_dict = {'T': t, 'V': v, 'sigma_v': sigma_v, 'sigma_t': sigma_t, 'delta_v': delta_v, 'delta_t': delta_t,
                  'ts_missing_val': ts_missing_val, 'tv_missing_val': tv_missing_val,
                  'aligned_read_count': aligned_read_count}

    # Second part

    data_ts_missing = pow(param_dict['delta_t'], param_dict['ts_missing_val'])
    data_tv_missing = pow(param_dict['delta_v'], param_dict['tv_missing_val'])

    print('calculating the likelihood nominator', file=sys.stderr)

    # todo could you double check that ll_nom is equivalent to this SQL statement:
    #  UPDATE ts_tv SET ll_nom = POWER(@delta_t, Ts) * POWER(@delta_v, Tv);

    init_ts_tv = pd.read_csv(ts_tv_file, names=['Taxon', 'Read_ID', 'Ts', 'Tv'], sep=',')
    init_ts_tv['ll_nom'] = init_ts_tv['Ts'].rpow(param_dict['delta_t']) * init_ts_tv['Tv'].rpow(param_dict['delta_v'])

    print(init_ts_tv, file=sys.stderr)

    total_taxa = pd.read_csv(taxa_file)
    total_taxa_count = len(total_taxa)

    print('calculating the proper likelihood', file=sys.stderr)

    init_ts_tv['Likelihood'] = np.nan

    # todo could you please double check that the following for loop is equivalent to this SQL statement:
    #  UPDATE ts_tv
    #  JOIN
    #  (SELECT read_id, sum(ll_nom)+(@total_taxa - count(taxa_id)) * @data_ts_missing * @data_tv_missing AS denominator
    #  FROM ts_tv GROUP BY read_id) AS den
    #  ON ts_tv.read_id = den.read_id
    #  SET Likelihood = ll_nom / denominator;

    for index, group in init_ts_tv.groupby('Read_ID'):
        init_ts_tv['Likelihood'] = init_ts_tv['ll_nom'].transform(
            lambda nom: nom / sum(nom) + ((total_taxa_count - len(group['Taxon'])) * data_ts_missing * data_tv_missing)
        )

    print(init_ts_tv, file=sys.stderr)

    # todo could you double check that the following dirichlet assignment is the equivalent to this SQL statement:
    #  UPDATE ts_tv SET Dirichlet_Assignment = CASE WHEN Likelihood >= 0.95 THEN 1 ELSE 0 END;

    # do the dirichlet assignment
    init_ts_tv['Dirichlet_Assignment'] = np.nan
    init_ts_tv['Dirichlet_Assignment'] = init_ts_tv.where(init_ts_tv['Likelihood'] >= 0.95, 0).astype('int')

    init_ts_tv.to_csv(output_matrix, sep=',', index=False)

    # TODO I found an encoder class so that the numpy integer types in the parameters dict could be
    #  converted to serializable types

    with open(output_params, 'w') as fout:
        json.dump(param_dict, fout, cls=NpEncoder)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    calculate_likelihoods(
        ts_tv_file=snakemake.input[0],
        readlen_file=snakemake.input[1],
        taxa_file=snakemake.input[2],
        config=snakemake.config,
        output_matrix=snakemake.output[0],
        output_params=snakemake.output[1]
    )
