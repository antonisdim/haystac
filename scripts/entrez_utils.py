#!/usr/bin/env python
# -*- coding: utf-8 -*-

import http.client
import sys
import time
import urllib.error
from datetime import datetime
import socket

from Bio import Entrez

# the maximum number of attempts to make for a failed query
MAX_RETRY_ATTEMPTS = 2

# time to wait in seconds before repeating a failed query
RETRY_WAIT_TIME = 2

ENTREZ_DB_NUCCORE = 'nuccore'
ENTREZ_DB_TAXA = 'taxonomy'

ENTREZ_RETMODE_XML = 'xml'
ENTREZ_RETMODE_TEXT = 'text'

ENTREZ_RETTYPE_FASTA = 'fasta'
ENTREZ_RETTYPE_GB = 'gb'

ENTREZ_RETMAX = 10 ** 9


def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def entrez_efetch(db, retmode, rettype, webenv, query_key, attempt=1):
    try:
        return Entrez.efetch(db=db,
                             retmode=retmode,
                             rettype=rettype,
                             retmax=ENTREZ_RETMAX,
                             webenv=webenv,
                             query_key=query_key)

    except http.client.HTTPException as e:
        print("Network problem: {}".format(e), file=sys.stderr)

        attempt += 1

        if attempt > MAX_RETRY_ATTEMPTS:
            print("Exceeded maximum attempts {}...".format(attempt), file=sys.stderr)
            return None
        else:
            time.sleep(RETRY_WAIT_TIME)
            print("Starting attempt {}...".format(attempt), file=sys.stderr)
            return entrez_efetch(db, retmode, rettype, webenv, query_key, attempt)

    except (http.client.IncompleteRead, urllib.error.URLError) as e:
        # TODO refactor this error handling -
        #  that should be self fixing as it an error related to the try except block of guts_of_entrez
        print("Ditching that batch", file=sys.stderr)
        print(e)
        return None


def guts_of_entrez(db, retmode, rettype, chunk, batch_size):
    # print info about number of records
    print("Downloading {} entries from NCBI {} database in batches of {} entries...\n"
          .format(len(chunk), db, batch_size), file=sys.stderr)
    # post NCBI query
    search_handle = Entrez.epost(db, id=",".join(map(str, chunk)))
    search_results = Entrez.read(search_handle)

    now = datetime.ctime(datetime.now())
    print("\t{} for a batch of {} records \n".format(now, len(chunk)), file=sys.stderr)

    handle = entrez_efetch(db, retmode, rettype, search_results["WebEnv"], search_results["QueryKey"])

    # print("got the handle", file=sys.stderr)
    if not handle:
        raise RuntimeError("The records from the following accessions could not be fetched: {}".format(','.join(chunk)))

    try:
        if retmode == ENTREZ_RETMODE_TEXT:
            yield handle.read()
        else:
            records = Entrez.read(handle)
            # print("got the records", file=sys.stderr)

            for rec in records:
                yield rec

    except (http.client.HTTPException, urllib.error.HTTPError, urllib.error.URLError,
            RuntimeError, Entrez.Parser.ValidationError, socket.error):

        for accession in chunk:
            try:
                guts_of_entrez(db, retmode, rettype, accession, batch_size)
            except (http.client.HTTPException, urllib.error.HTTPError, urllib.error.URLError,
                    RuntimeError, Entrez.Parser.ValidationError, socket.error):
                print("Ditching this accession as it is a bad record {}.".format(accession), file=sys.stderr)

        # TODO refactor this error handling
