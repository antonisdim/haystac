#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import sys
import pandas as pd
import itertools
import time
import urllib.error

from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import guts_of_entrez, ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_GB

# todo should it in the config ? needs to be the same as in entrez.smk file. Feel like I can do it better
CHUNK_SIZE = 20
TOO_MANY_REQUESTS_WAIT = 7
MAX_RETRY_ATTEMPTS = 10


def gen_dict_extract(key, var):
    """Find keys in nested dictionaries"""
    if hasattr(var, 'items'):
        for k, v in var.items():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result
    elif isinstance(var, list):
        for d in var:
            for result in gen_dict_extract(key, d):
                yield result


def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def entrez_nuccore_query(input_file, config, query_chunk_num, output_file, attempt=1):
    """
    Query the NCBI nuccore database to get a list of sequence accessions and their metadata.
    """

    time.sleep(int(query_chunk_num) // 3)
    Entrez.email = config['entrez_email']

    accessions = pd.read_csv(input_file, sep='\t').squeeze().to_list()
    entrez_query_list = next(itertools.islice(chunker(accessions, CHUNK_SIZE), int(query_chunk_num), None))
    entrez_query_list = [acc + '[Accession]' for acc in entrez_query_list]

    entrez_query = ' OR '.join(entrez_query_list)

    # entrez_query = config['entrez']['queries'][query_chunk]

    # get list of entries for given query
    print("Getting list of Accessions for term={} ...\n".format(entrez_query), file=sys.stderr)

    retmax = config['entrez_batchsize']
    counter = 0
    dictwriter_counter = 0

    with open(output_file, 'a') as fout:
        fileEmpty = os.stat(output_file).st_size == 0
        fieldnames = ['GBSeq_accession-version', 'TSeq_taxid', 'GBSeq_organism', 'GBSeq_definition', 'GBSeq_length']
        w = csv.DictWriter(fout, fieldnames, delimiter='\t', extrasaction="ignore")
        if fileEmpty:
            w.writeheader()

        while True:
            handle = Entrez.esearch(db=ENTREZ_DB_NUCCORE, term=entrez_query, retmax=retmax, idtype="acc",
                usehistory='y', retstart=retmax * counter, rettype=ENTREZ_RETTYPE_GB,
                retmode=ENTREZ_RETMODE_XML)

            handle_reader = Entrez.read(handle)
            accessions = handle_reader['IdList']

            if counter == 0:
                resultset = int(handle_reader['Count'])
                total_records = resultset
                print('Total number of sequences {}\n'.format(resultset), file=sys.stderr)
            else:
                print('Remaining sequences to be downloaded {}\n'.format(resultset), file=sys.stderr)

            if not accessions:
                # stop iterating when we get an empty resultset
                if dictwriter_counter == total_records:
                    print("A total of {} records have been saved successfully.\n".format(total_records),
                        file=sys.stderr)
                else:
                    print('dictwriter_counter\t', dictwriter_counter)
                    print('total_records\t', total_records)
                    raise RuntimeError("A total of {} records have been saved successfully. Please check the relevant "
                                       "log file to see which ones failed.\n".format(dictwriter_counter))
                break

            try:
                records = guts_of_entrez(ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_GB, accessions,
                    retmax)
                for node in records:
                    # print(node)

                    node['TSeq_taxid'] = [val.replace('taxon:', '')
                                          for val in gen_dict_extract('GBQualifier_value', node) if
                                          'taxon:' in val].pop()
                    # print("iterating on node", file=sys.stderr)
                    print(node)
                    w.writerow(node)
                    dictwriter_counter += 1

                print("done for this slice\n", file=sys.stderr)

                counter += 1
                resultset -= retmax
                if resultset < 0:
                    resultset = 0

            except urllib.error.URLError as e:
                print("Network problem: {}".format(e), file=sys.stderr)

                attempt += 1

                if attempt > MAX_RETRY_ATTEMPTS:
                    print("Exceeded maximum attempts {}...".format(attempt), file=sys.stderr)
                    return None
                else:
                    time.sleep(TOO_MANY_REQUESTS_WAIT)
                    print("Starting attempt {}...".format(attempt), file=sys.stderr)
                    entrez_nuccore_query(input_file, config, query_chunk_num, output_file, attempt=1)

        print("COMPLETE", file=sys.stderr)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_nuccore_query(
        input_file=snakemake.input[0],
        config=snakemake.config,
        query_chunk_num=snakemake.wildcards.chunk,
        output_file=snakemake.output[0]
    )
