#!/usr/bin/env python
# -*- coding: utf-8 -*-
from datetime import datetime
import sys
import time

from Bio import Entrez

import http.client
import urllib.error

from socket import error as socketerror


def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def get_entrez_query(config, query, output_file):
    """
    Query the NCBI Entraz database to get a list of sequence accessions and their taxa details for the given query.
    """

    Entrez.email = config['entrez']['email']

    # get list of entries for given query
    print("Getting list of GIs for term={} ...\n".format(query))

    handle = Entrez.esearch(db=config['entrez']['db'], term=query, retmax=config['entrez']['retmax'], idtype="acc")
    gi_list = Entrez.read(handle)['IdList']

    for acc_num in chunker(gi_list, 100):
        print(acc_num)

        # print info about number of proteins
        if config['verbose_output']:
            print("Downloading %s entries from NCBI %s database in batches of %s entries...\n"
                  .format(len(acc_num), config['entrez']['db'], config['entrez']['batchSize']))

        # post NCBI query
        search_handle = Entrez.epost(config['entrez']['db'], id=",".join(acc_num))
        search_results = Entrez.read(search_handle)
        webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]

        for start in range(0, len(acc_num), config['entrez']['batchSize']):
            # print info
            tnow = datetime.now()
            print("\t%s\t%s / %s\n" % (datetime.ctime(tnow), start, len(acc_num)))
            try:
                handle = Entrez.efetch(db=config['entrez']['db'],
                                       retmode=config['entrez']['retmode'],
                                       rettype=config['entrez']['rettype'],
                                       retstart=start,
                                       retmax=config['entrez']['batchSize'],
                                       webenv=webenv,
                                       query_key=query_key)

            except http.client.HTTPException as e:
                print("Network problem: {}".format(e))
                print("Second (and final) attempt...")

                handle = Entrez.efetch(db=config['entrez']['db'],
                                       retmode=config['entrez']['retmode'],
                                       rettype=config['entrez']['rettype'],
                                       retstart=start,
                                       retmax=config['entrez']['batchSize'],
                                       webenv=webenv,
                                       query_key=query_key)

            except (http.client.IncompleteRead, urllib.error.URLError):
                print("Ditching that batch")
                continue

            print("got the handle")

            try:
                # records = Entrez.read(handle)
                records = Entrez.read(handle)
                print("got the records")

            except (http.client.HTTPException, urllib.error.HTTPError, urllib.error.URLError,
                    RuntimeError, Entrez.Parser.ValidationError, socketerror) as e:
                # print "Network problem: %s" % e
                # print "Second (and final) attempt..."
                # records = Entrez.read(handle)
                # print "got the handle"
                print("Ditching that batch")
                continue

            # print records

            for node in records:

                print("iterating on node")

                record = dict()
                record["code"] = node["TSeq_accver"]
                record["taxid"] = node["TSeq_taxid"]
                record["length"] = node["TSeq_length"]
                record["org_name"] = node["TSeq_orgname"]
                # record["sequence"] = node["TSeq_sequence"]

            time.sleep(10)

        print("done for this slice")


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    get_entrez_query(
        config=snakemake.config,
        query=snakemake.wildcard.query,
        output_file=snakemake.output[0]
    )
