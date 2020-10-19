#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


from xml.etree import ElementTree

import http.client
import requests
import socket
import sys
import time
import urllib.error
from Bio import Entrez
from datetime import datetime
from urllib.parse import quote_plus

# base url of the Entrez web service
ENTREZ_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# maximum 1 request per X seconds
ENTREZ_RATE_LIMIT = 1

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


# TODO refactor this out
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


# TODO refactor this out
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


# TODO refactor this out
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


def entrez_request(url):
    """
    Helper function to ensure that we never exceed the rate limit.
    """
    # TODO if config['api_key'] is set then append it here...
    r = requests.get(ENTREZ_URL + url)

    if not r.ok:
        r.raise_for_status()

    time.sleep(ENTREZ_RATE_LIMIT)

    return r


def entrez_esearch(database, query):
    """
    Execute an Entrez esearch query and return the search keys
    """
    r = entrez_request(f"esearch.fcgi?db={database}&term={quote_plus(query)}&usehistory=y")

    # parse the XML result
    etree = ElementTree.XML(r.text)

    # get the search keys
    key = etree.find("QueryKey").text
    webenv = etree.find("WebEnv").text
    id_list = [id.text for id in etree.findall(".//Id")]

    return key, webenv, id_list


def entrez_esummary(database, key, webenv):
    """
    Fetch the Entrez esummary records for an esearch query.
    """
    r = entrez_request(f"esummary.fcgi?db={database}&query_key={key}&WebEnv={webenv}")

    return ElementTree.XML(r.text)


def entrez_range_accessions(accession, first, last):
    """
    Return a range between two accession codes (e.g. ABC001 -> ABC005)
    """
    # sanity check that the items are equal length and in order
    assert len(first) == len(last) and first < last

    # find the index of the first difference
    idx = [i for i in range(len(first)) if first[i] != last[i]][0]

    try:
        # return the range
        return [f"{first[:idx]}{item}" for item in range(int(first[idx:]), int(last[idx:]) + 1)]
    except ValueError:
        print(
            f"ERROR: Could not resolve the accession range {first}-{last} for master record {accession}",
            file=sys.stderr,
        )
        exit(1)
