#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import http.client
import sys
import urllib.error
from datetime import datetime
from socket import error as socketerror

from Bio import Entrez
from entrez_utils import chunker, entrez_efetch


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
            print(acc_num, file=sys.stderr)

            # print info about number of proteins
            print("Downloading {} entries from NCBI {} database in batches of {} entries...\n"
                  .format(len(acc_num), 'nuccore', config['entrez']['batchSize']), file=sys.stderr)

            # post NCBI query
            search_handle = Entrez.epost('nuccore', id=",".join(acc_num))
            search_results = Entrez.read(search_handle)

            for start in range(0, len(acc_num), config['entrez']['batchSize']):
                # print info
                tnow = datetime.now()
                print("\t{}\t{} / {}\n".format(datetime.ctime(tnow), start, len(acc_num)), file=sys.stderr)

                handle = entrez_efetch(config, 'nuccore', start, search_results["WebEnv"], search_results["QueryKey"])

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
