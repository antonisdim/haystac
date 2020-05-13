#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import pandas as pd

import os.path


def entrez_pick_sequences(config, nuccore_file, taxa_file, output_file, query):
    accessions = pd.read_csv(nuccore_file, sep='\t')
    taxa = pd.read_csv(taxa_file, sep='\t')
    rank = config['entrez_rank']

    print("read the accessions and the taxa", file=sys.stderr)

    sequences = pd.merge(accessions, taxa, on=['TSeq_taxid'], how='outer')

    sequences = sequences[~sequences[rank].isnull()]

    selected_sequences = sequences.loc[
        sequences.groupby(rank)['GBSeq_length'].idxmax(), ['species', 'GBSeq_accession-version']]

    print(selected_sequences, file=sys.stderr)

    if os.path.isfile("database_inputs/prok_representative_genomes.txt"):
        refseq_genomes = pd.read_csv("{query}/entrez/{query}-refseq-genomes.tsv".format(query=query), sep='\t')
        genbank_genomes = pd.read_csv("{query}/entrez/{query}-genbank-genomes.tsv".format(query=query), sep='\t')
        assemblies = pd.read_csv("{query}/entrez/{query}-assemblies.tsv".format(query=query), sep='\t')
        refseq_plasmids = pd.read_csv("{query}/entrez/{query}-refseq-plasmids.tsv".format(query=query), sep='\t')
        genbank_plasmids = pd.read_csv("{query}/entrez/{query}-genbank-plasmids.tsv".format(query=query), sep='\t')

        # the entrez query might give a different accession for a certain species than the refseq rep one and
        # I don't want that. If the species exists in the refseq I want to keep the refseq records
        selected_sequences = selected_sequences[~selected_sequences.species.isin(refseq_genomes.species)]
        selected_sequences = selected_sequences[~selected_sequences.species.isin(genbank_genomes.species)]
        selected_sequences = selected_sequences[~selected_sequences.species.isin(assemblies.species)]
        selected_sequences = selected_sequences[~selected_sequences.species.isin(refseq_plasmids.species)]
        selected_sequences = selected_sequences[~selected_sequences.species.isin(genbank_plasmids.species)]

    print("selected the longest sequence per species, writing it to a file", file=sys.stderr)

    selected_sequences.to_csv(output_file, sep='\t', header=True, index=False)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_pick_sequences(
        config=snakemake.config,
        nuccore_file=snakemake.input[0],
        taxa_file=snakemake.input[1],
        output_file=snakemake.output[0],
        query=snakemake.wildcards.query,
    )
