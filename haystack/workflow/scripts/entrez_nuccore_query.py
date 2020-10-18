#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import csv

import sys
import json
import time
import pandas as pd

import requests
from urllib.parse import quote_plus
from xml.etree import cElementTree as ElementTree

# base url of the Entrez web service
ENTREZ_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# number of seconds to wait before making a new request
ENTREZ_WAIT_TIME = 1

# how many seconds to wait between retries (weighted by number of attempts)
ENTREZ_RETRY_WAIT = 2

# maximum number of times to retry fetching an Entrez record before giving up
ENTREZ_MAX_RETRY = 3


def entrez_esearch(database, query, attempts=1):
    """
    Execute an Entrez esearch query and return the search keys
    """
    r = requests.get(ENTREZ_URL + f"esearch.fcgi?db={database}&term={quote_plus(query)}&usehistory=y")

    if not r.ok:
        # handle Too Many Requests error
        if attempts < ENTREZ_MAX_RETRY:
            wait = ENTREZ_RETRY_WAIT * attempts
            print(
                "WARNING: Entrez esearch query failed on attempt #{attempt}, retrying after {wait} seconds.",
                file=sys.stderr,
            )
            time.sleep(wait)
            return entrez_esearch(database, query, attempts + 1)
        else:
            r.raise_for_status()

    etree = ElementTree.XML(r.text)

    # get the search keys
    key = etree.find("QueryKey").text
    webenv = etree.find("WebEnv").text

    return key, webenv


def entrez_esummary(database, key, webenv, attempts=1):
    """
    Fetch the Entrez esummary records for an esearch query.
    """
    r = requests.get(ENTREZ_URL + f"esummary.fcgi?db={database}&query_key={key}&WebEnv={webenv}")

    if not r.ok:
        # handle Too Many Requests error
        if attempts < ENTREZ_MAX_RETRY:
            wait = ENTREZ_RETRY_WAIT * attempts
            print(
                "WARNING: Entrez esummary query failed on attempt #{attempt}, retrying after {wait} seconds.",
                file=sys.stderr,
            )
            time.sleep(wait)
            return entrez_esummary(database, key, webenv, attempts + 1)
        else:
            r.raise_for_status()

    return ElementTree.XML(r.text)


def element_tree_to_dict(etree):
    """
    Convert an ElementTree object into a list of dicts.
    """
    data = []

    for row_node in etree:
        row = {}
        for col_node in row_node:
            col = col_node.attrib.get("Name", col_node.tag)
            row[col] = col_node.text

        data.append(row)

    return data


def entrez_nuccore_query(query, output_file):
    """
    Query the NCBI nuccore database to get a list of sequence accessions and their metadata.
    """
    # execute the search
    key, webenv = entrez_esearch("nuccore", query)

    # fetch the results
    etree = entrez_esummary("nuccore", key, webenv)

    # convert the ElementTree into a a list of dicts
    data = element_tree_to_dict(etree)

    # print(data)

    with open(output_file, "w") as fout:
        w = csv.DictWriter(fout, data[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(data)


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    entrez_nuccore_query(
        query=snakemake.config["query"], output_file=snakemake.output[0],
    )
