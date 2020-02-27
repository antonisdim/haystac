#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import sys

from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import guts_of_entrez, ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_FASTA


def entrez_nuccore_query(config, query, output_file):
    """
    Query the NCBI nuccore database to get a list of sequence accessions and their metadata.
    """
    Entrez.email = config['entrez']['email']
    entrez_query = config['entrez']['queries'][query]

    # get list of entries for given query
    print("Getting list of GIs for term={} ...\n".format(entrez_query), file=sys.stderr)

    retmax = config['entrez']['batchSize']
    counter = 0

    while True:
        handle = Entrez.esearch(db=ENTREZ_DB_NUCCORE, term=entrez_query, retmax=retmax, idtype="acc",
                                usehistory='y', retstart=retmax*counter)
        accessions = Entrez.read(handle)['IdList']

        if not accessions:
            # stop iterating when we get an empty resultset
            break

        with open(output_file, 'w') as fout:
            fieldnames = ['TSeq_accver', 'TSeq_taxid', 'TSeq_orgname', 'TSeq_defline', 'TSeq_length']
            w = csv.DictWriter(fout, fieldnames, delimiter='\t', extrasaction="ignore")
            w.writeheader()

            # TODO can we not download the FASTA sequence, as this can be many MB for each result
            records = guts_of_entrez(ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_FASTA,
                                     accessions, retmax)

            for node in records:
                print("iterating on node", file=sys.stderr)
                w.writerow(node)

            print("done for this slice", file=sys.stderr)

        counter += 1

        print("COMPLETE", file=sys.stderr)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_nuccore_query(
        config=snakemake.config,
        query=snakemake.wildcards.query,
        output_file=snakemake.output[0]
    )
