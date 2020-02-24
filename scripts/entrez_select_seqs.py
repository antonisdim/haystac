#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import pandas as pd


def entrez_select_taxa_seqs(nuccore_file, taxa_file, output_file):

    accessions = pd.read_csv(nuccore_file, sep='\t')
    taxa = pd.read_csv(taxa_file, sep='\t')

    print("read the accessions and the taxa", file=sys.stderr)

    # TODO make rank configurable in the config.yaml file
    rank = 'species'

    taxaccessions = pd.merge(accessions, taxa, on=['TSeq_taxid'], how='outer')

    # TODO if rank==subspecies, you must issue a warning that the resultset sizes are different
    taxaccessions = taxaccessions[~taxaccessions[rank].isnull()]

    # TODO groupby on species name but species are not unique
    selected_sequences = taxaccessions.loc[
        taxaccessions.groupby(rank)['TSeq_length'].idxmax(), ['TSeq_orgname', 'TSeq_accver']]

    print("selected the longest sequence per species, writing it to a file", file=sys.stderr)

    selected_sequences.to_csv(output_file, sep='\t', header=True, index=False)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_select_taxa_seqs(
        nuccore_file=snakemake.input[0],
        taxa_file=snakemake.input[1],
        output_file=snakemake.output[0]
    )
