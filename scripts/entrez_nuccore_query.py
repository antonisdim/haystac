#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import sys

from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import chunker, guts_of_entrez


def entrez_nuccore_query(config, query, output_file):
    """
    Query the NCBI nuccore database to get a list of sequence accessions and their metadata.
    """
    Entrez.email = config['entrez']['email']
    entrez_query = config['entrez']['queries'][query]

    # get list of entries for given query
    print("Getting list of GIs for term={} ...\n".format(entrez_query), file=sys.stderr)

    handle = Entrez.esearch(db='nuccore', term=entrez_query, retmax=config['entrez']['retmax'], idtype="acc")
    accessions = Entrez.read(handle)['IdList']

    with open(output_file, 'w') as fout:
        fieldnames = ['TSeq_accver', 'TSeq_taxid', 'TSeq_orgname', 'TSeq_defline', 'TSeq_length']
        w = csv.DictWriter(fout, fieldnames, delimiter='\t', extrasaction="ignore")
        w.writeheader()

        for acc_num in chunker(accessions, 100):

            records = guts_of_entrez('nuccore', acc_num, config)

            for node in records:
                print("iterating on node", file=sys.stderr)
                w.writerow(node)

            print("done for this slice", file=sys.stderr)

    print("COMPLETE", file=sys.stderr)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_nuccore_query(
        config=snakemake.config,
        query=snakemake.wildcards.query,
        output_file=snakemake.output[0]
    )
