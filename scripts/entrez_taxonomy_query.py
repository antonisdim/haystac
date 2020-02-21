#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import sys

import pandas as pd
from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import chunker, guts_of_entrez


def entrez_taxonomy_query(config, nuccore_file, output_file):
    """
    Query the NCBI taxonomy database to get a taxa details for all nuccore sequences.
    """

    Entrez.email = config['entrez']['email']

    # load the unique list of taxa from the nuccore resultset
    accessions = pd.read_csv(nuccore_file, sep='\t', usecols=['TSeq_taxid'], squeeze=True).unique()

    with open(output_file, 'w') as fout:
        fieldnames = ['TaxId', 'genus', 'family', 'species', 'subspecies']
        w = csv.DictWriter(fout, fieldnames, delimiter='\t', extrasaction="ignore")
        w.writeheader()

        for chunk in chunker(accessions, 100):

            records = guts_of_entrez('taxonomy', chunk, config)

            for node in records:
                taxon = dict()
                taxon['TaxId'] = node['TaxId']
                for item in node['LineageEx']:
                    if item["Rank"] in fieldnames:
                        taxon[item["Rank"]] = item['ScientificName']

                if 'species' not in taxon.keys():
                    if node['Rank'] == 'species':
                        taxon['species'] = node['ScientificName']

                # print(taxon)

                w.writerow(taxon)

                print("done for this slice", file=sys.stderr)

    print("COMPLETE", file=sys.stderr)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_taxonomy_query(
        config=snakemake.config,
        nuccore_file=snakemake.input[0],
        output_file=snakemake.output[0]
    )
