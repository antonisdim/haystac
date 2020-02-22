#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd


def entrez_select_taxa_seqs(accesions_file, taxa_file, output_file):
    accessions = pd.read_csv(accesions_file, sep='\t', skiprows=[0], header=None,
                             names=['accession', 'taxid', 'org_name', 'defline', 'length'], dtype={'taxid': 'str'})

    taxa = pd.read_csv(taxa_file, sep='\t', skiprows=[0],
                       header=None, names=['taxid', 'genus', 'family', 'species', 'subspecies'], dtype={'taxid': 'str'})

    print("read the accessions and the taxa", file=sys.stderr)

    taxaccessions = pd.merge(accessions, taxa, on=['taxid'], how='outer')

    taxaccessions = taxaccessions[~taxaccessions['species'].isnull()]

    selected_sequences = taxaccessions.loc[
        taxaccessions.groupby('species')['length'].idxmax(), ['species', 'accession']]

    print("selected the longest sequence per species, writing it to a file", file=sys.stderr)

    selected_sequences.to_csv(output_file, sep='\t', header=False, index=False)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_select_taxa_seqs(
        accesions_file=snakemake.input[0],
        taxa_file=snakemake.input[1],
        output_file=snakemake.output[0]
    )
