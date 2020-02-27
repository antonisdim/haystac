#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import pandas as pd


def entrez_pick_sequences(config, nuccore_file, taxa_file, output_file):

    accessions = pd.read_csv(nuccore_file, sep='\t')
    taxa = pd.read_csv(taxa_file, sep='\t')

    print("read the accessions and the taxa", file=sys.stderr)

    taxaccessions = pd.merge(accessions, taxa, on=['TSeq_taxid'], how='outer')

    taxaccessions = taxaccessions[~taxaccessions[config['rank']].isnull()]

    selected_sequences = taxaccessions.loc[
        taxaccessions.groupby(config['rank'])['TSeq_length'].idxmax(), ['TSeq_orgname', 'TSeq_accver']]

    print("selected the longest sequence per species, writing it to a file", file=sys.stderr)

    selected_sequences.to_csv(output_file, sep='\t', header=True, index=False)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_pick_sequences(
        config=snakemake.config,
        nuccore_file=snakemake.input[0],
        taxa_file=snakemake.input[1],
        output_file=snakemake.output[0]
    )
