#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import http.client
import socket
import sys
import time
import urllib.error
from Bio import Entrez
from datetime import datetime

# the maximum number of attempts to make for a failed query
MAX_RETRY_ATTEMPTS = 2

# time to wait in seconds before repeating a failed query
RETRY_WAIT_TIME = 2

TOO_MANY_REQUESTS_WAIT = 20

ENTREZ_DB_NUCCORE = "nuccore"
ENTREZ_DB_TAXA = "taxonomy"
ENTREZ_DB_ASSEMBLY = "assembly"

ENTREZ_RETMODE_XML = "xml"
ENTREZ_RETMODE_TEXT = "text"

ENTREZ_RETTYPE_FASTA = "fasta"
ENTREZ_RETTYPE_GB = "gb"

ENTREZ_RETMAX = 10 ** 9


def chunker(seq, size):
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def entrez_efetch(db, retmode, rettype, webenv, query_key, attempt=1):
    try:
        return Entrez.efetch(
            db=db, retmode=retmode, rettype=rettype, retmax=ENTREZ_RETMAX, webenv=webenv, query_key=query_key,
        )

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
        attempt += 1

        if attempt > MAX_RETRY_ATTEMPTS:
            print("Exceeded maximum attempts {}...".format(attempt), file=sys.stderr)
            return None
        elif e.code == 429:
            time.sleep(RETRY_WAIT_TIME)
            print("Starting attempt {}...".format(attempt), file=sys.stderr)
            return entrez_efetch(db, retmode, rettype, webenv, query_key, attempt)
        else:
            print("Discarding that batch", file=sys.stderr)
            print(e, file=sys.stderr)
            return None


def guts_of_entrez(db, retmode, rettype, chunk, batch_size):
    # print info about number of records
    print(
        "Downloading {} entries from the NCBI {} database in batches of {} entries...\n".format(
            len(chunk), db, batch_size
        ),
        file=sys.stderr,
    )
    # post NCBI query
    search_handle = Entrez.epost(db, id=",".join(map(str, chunk)))
    search_results = Entrez.read(search_handle)

    now = datetime.ctime(datetime.now())
    print("\t{} for a batch of {} records \n".format(now, len(chunk)), file=sys.stderr)

    handle = entrez_efetch(db, retmode, rettype, search_results["WebEnv"], search_results["QueryKey"])

    # print("got the handle", file=sys.stderr)
    if not handle:
        raise RuntimeError("The records from the following accessions could not be fetched: {}".format(",".join(chunk)))

    try:
        if retmode == ENTREZ_RETMODE_TEXT:
            yield handle.read()
        else:
            records = Entrez.read(handle)
            # print("got the records", file=sys.stderr)

            for rec in records:
                yield rec

    except (
        http.client.HTTPException,
        urllib.error.HTTPError,
        urllib.error.URLError,
        RuntimeError,
        Entrez.Parser.ValidationError,
        socket.error,
    ):

        for accession in chunk:
            try:
                yield guts_of_entrez(db, retmode, rettype, accession, batch_size)
            except (
                http.client.HTTPException,
                urllib.error.HTTPError,
                urllib.error.URLError,
                RuntimeError,
                Entrez.Parser.ValidationError,
                socket.error,
            ):
                print(
                    "Discarding this accession as it is a bad record {}.".format(accession), file=sys.stderr,
                )


def get_accession_ftp_path(accession, config, attempt=1):
    """Get a valid NCBI ftp path from an accession."""

    Entrez.email = config["email"]
    try:
        handle = Entrez.esearch(db=ENTREZ_DB_ASSEMBLY, term=accession + ' AND "latest refseq"[filter]')
        # or handle = Entrez.esearch(db=ENTREZ_DB_ASSEMBLY,
        # term=accession + ' AND ((latest[filter] OR "latest refseq"[filter])')
        assembly_record = Entrez.read(handle)
        esummary_handle = Entrez.esummary(db=ENTREZ_DB_ASSEMBLY, id=assembly_record["IdList"], report="full")
        try:
            esummary_record = Entrez.read(esummary_handle, validate=False)
        except RuntimeError:
            return ""
        refseq_ftp = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
        genbank_ftp = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_GenBank"]

        if refseq_ftp != "":
            return refseq_ftp
        else:
            return genbank_ftp

    except urllib.error.HTTPError as e:
        if e.code == 429:

            attempt += 1

            if attempt > MAX_RETRY_ATTEMPTS:
                print("Exceeded maximum attempts {}...".format(attempt), file=sys.stderr)
                return None
            else:
                time.sleep(TOO_MANY_REQUESTS_WAIT)
                get_accession_ftp_path(accession, config, attempt)

        else:
            raise RuntimeError("There was a urllib.error.HTTPError with code {}".format(e))

    except IndexError:
        time.sleep(TOO_MANY_REQUESTS_WAIT)
        get_accession_ftp_path(accession, config, attempt)


