#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import sys

import pandas as pd
from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import chunker, guts_of_entrez, ENTREZ_DB_TAXA, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_FASTA


def entrez_taxonomy_query(config, nuccore_file, output_file):
    """
    Query the NCBI taxonomy database to get a taxa details for all nuccore sequences.
    """

    Entrez.email = config['entrez']['email']

    # load the unique list of taxa from the nuccore resultset
    accessions = pd.read_csv(nuccore_file, sep='\t', usecols=['TSeq_taxid'], squeeze=True).unique()

    with open(output_file, 'w') as fout:
        columns = ['TSeq_taxid', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
        w = csv.DictWriter(fout, columns, delimiter='\t', extrasaction="ignore")
        w.writeheader()

        resultset = len(accessions)

        for chunk in chunker(accessions, config['entrez']['batchSize']):
            print('Remaining sequences to have their taxids downloaded {}\n'.format(resultset), file=sys.stderr)
            records = guts_of_entrez(ENTREZ_DB_TAXA, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_FASTA,
                                     chunk, config['entrez']['batchSize'])

            for node in records:
                taxon = dict()
                taxon['TSeq_taxid'] = node['TaxId']
                for item in node['LineageEx']:
                    if item["Rank"] in columns:
                        taxon[item["Rank"]] = item['ScientificName']

                if node['Rank'] in columns:
                    taxon[node['Rank']] = node['ScientificName']

                if taxon['species'] and not taxon.get('subspecies'):
                    taxon['subspecies'] = taxon['species'] + ' ssp.'

                w.writerow(taxon)

            resultset -= config['entrez']['batchSize']
            print("done for this slice\n", file=sys.stderr)

    print("COMPLETE", file=sys.stderr)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_taxonomy_query(
        config=snakemake.config,
        nuccore_file=snakemake.input[0],
        output_file=snakemake.output[0]
    )
