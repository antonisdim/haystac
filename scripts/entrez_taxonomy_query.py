#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import http.client
import sys
import urllib.error
from datetime import datetime
from socket import error as socketerror

import pandas as pd
from Bio import Entrez
from entrez_utils import chunker, entrez_efetch


def entrez_taxonomy_query(config, nuccore_file, output_file):
    """
    Query the NCBI taxonomy database to get a taxa details for all nuccore sequences.
    """

    Entrez.email = config['entrez']['email']

    # load the unique list of taxa from the nuccore resultset
    accessions = pd.read_csv(nuccore_file, sep='\t', usecols=['TSeq_taxid'], squeeze=True).unique()

    with open(output_file, 'w') as fout:
        fieldnames = ['genus', 'family', 'species', 'subspecies']
        w = csv.DictWriter(fout, fieldnames, delimiter='\t', extrasaction="ignore")
        w.writeheader()

        for chunk in chunker(accessions, 100):
            # print info about number of records
            print("Downloading {} entries from NCBI {} database in batches of {} entries...\n"
                  .format(len(accessions), 'taxonomy', config['entrez']['batchSize']), file=sys.stderr)

            # post NCBI query
            search_handle = Entrez.epost('taxonomy', id=",".join(map(str, chunk)))
            search_results = Entrez.read(search_handle)

            for start in range(0, len(chunk), config['entrez']['batchSize']):
                # print info
                tnow = datetime.now()
                print("\t{}\t{} / {}\n".format(datetime.ctime(tnow), start, len(chunk)), file=sys.stderr)

                handle = entrez_efetch(config, 'taxonomy', start, search_results["WebEnv"], search_results["QueryKey"])

                if not handle:
                    continue

                print("got the handle", file=sys.stderr)

                try:
                    # records = Entrez.read(handle)
                    records = Entrez.read(handle)
                    print("got the records", file=sys.stderr)

                except (http.client.HTTPException, urllib.error.HTTPError, urllib.error.URLError,
                        RuntimeError, Entrez.Parser.ValidationError, socketerror) as e:
                    print("Ditching that batch", file=sys.stderr)
                    continue

                print(records)

                for node in records:
                    taxon = dict()
                    taxon['TaxId'] = node['TaxId']
                    for item in node['LineageEx']:
                        if item["Rank"] in fieldnames:
                            taxon[item["Rank"]] = item['ScientificName']

                    if 'species' not in taxon.keys():
                        if node['Rank'] == 'species':
                            taxon['species'] = node['ScientificName']

                    print(taxon)

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
