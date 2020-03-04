#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import sys

from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import guts_of_entrez, ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_GB


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


def entrez_nuccore_query(config, query, output_file):
    """
    Query the NCBI nuccore database to get a list of sequence accessions and their metadata.
    """
    Entrez.email = config['entrez']['email']
    entrez_query = config['entrez']['queries'][query]

    # get list of entries for given query
    print("Getting list of Accessions for term={} ...\n".format(entrez_query), file=sys.stderr)

    retmax = config['entrez']['batchSize']
    counter = 0
    dictwriter_counter = 0

    with open(output_file, 'a') as fout:
        fieldnames = ['GBSeq_accession-version', 'TSeq_taxid', 'GBSeq_organism', 'GBSeq_definition', 'GBSeq_length']
        w = csv.DictWriter(fout, fieldnames, delimiter='\t', extrasaction="ignore")
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

            print('Remaining sequences to be downloaded {}\n'.format(resultset), file=sys.stderr)

            if not accessions:
                # stop iterating when we get an empty resultset
                if dictwriter_counter == total_records:
                    print("A total of {} records have been saved successfully.\n".format(total_records),
                          file=sys.stderr)
                else:
                    print("WARNING: A total of {} records have been saved successfully. "
                          "Please check the relevant log file to see which ones failed.\n".format(total_records),
                          file=sys.stderr)
                break

            records = guts_of_entrez(ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_GB,
                                     accessions, retmax)

            for node in records:
                # print(node)

                node['TSeq_taxid'] = [val.replace('taxon:', '')
                                      for val in gen_dict_extract('GBQualifier_value', node) if
                                      'taxon:' in val].pop()
                # print("iterating on node", file=sys.stderr)

                w.writerow(node)
                dictwriter_counter += 1

            print("done for this slice\n", file=sys.stderr)

            counter += 1
            resultset -= retmax
            if resultset < 0:
                resultset = 0

        print("COMPLETE", file=sys.stderr)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_nuccore_query(
        config=snakemake.config,
        query=snakemake.wildcards.query,
        output_file=snakemake.output[0]
    )
