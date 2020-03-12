#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sys


def calculate_prob_model_params(ts_tv_file, readlen_file, config):
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

    parameters = {'T': t, 'V': v, 'sigma_v': sigma_v, 'sigma_t': sigma_t, 'delta_v': delta_v, 'delta_t': delta_t,
                  'ts_missing_val': ts_missing_val, 'tv_missing_val': tv_missing_val,
                  'aligned_read_count': aligned_read_count}

    return parameters


def calculate_likelihoods(ts_tv_file, readlen_file, taxa_file, config, output_matrix, output_params):
    param_dict = calculate_prob_model_params(ts_tv_file, readlen_file, config)

    data_ts_missing = pow(param_dict['delta_t'], param_dict['ts_missing_val'])
    data_tv_missing = pow(param_dict['delta_v'], param_dict['tv_missing_val'])

    print('calculating the likelihood nominator', file=sys.stderr)

    init_ts_tv = pd.read_csv(ts_tv_file, names=['Taxon', 'Read_ID', 'Ts', 'Tv'], sep=',')
    init_ts_tv['ll_nom'] = init_ts_tv['Ts'].rpow(param_dict['delta_t']) * init_ts_tv['Tv'].rpow(param_dict['delta_v'])

    print(init_ts_tv, file=sys.stderr)

    total_taxa = pd.read_csv(taxa_file, sep=',')
    total_taxa_count = len(total_taxa)

    print('calculating the proper likelihood', file=sys.stderr)

    init_ts_tv['Likelihood'] = np.nan

    for index, group in init_ts_tv.groupby('Read_ID'):
        init_ts_tv['Likelihood'] = init_ts_tv['ll_nom'].transform(lambda x: x / sum(x) +
                                ((total_taxa_count - len(group['Taxon'])) * data_ts_missing * data_tv_missing))

    print(init_ts_tv, file=sys.stderr)

    # do the dirichlet assignment

    init_ts_tv['Dirichlet_Assignment'] = np.nan
    init_ts_tv['Dirichlet_Assignment'] = init_ts_tv.where(init_ts_tv['Likelihood'] >= 0.95, 0).astype('int')

    init_ts_tv.to_csv(output_matrix, sep=',', index=False)

    pd.DataFrame.from_dict(param_dict, orient='index', columns=['Parameters']).to_csv(output_params, sep=',')


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    calculate_likelihoods(
        ts_tv_file=snakemake.input[0],
        readlen_file=snakemake.input[1],
        taxa_file=snakemake.input[2],
        config=snakemake.config,
        output_matrix=snakemake.output[0],
        output_params=snakemake.output[1])
