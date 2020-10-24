#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from xml.etree import ElementTree

import os
import requests
import sys
import time
import yaml
from urllib.parse import urlencode

# base url of the Entrez web service
ENTREZ_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# tool registration with NCBI, see https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Frequency_Timing_and_Registrati
ENTREZ_TOOL = "haystack"
ENTREZ_EMAIL = "antonisdim41@gmail.com"

# maximum 1 request per X seconds
ENTREZ_WAIT_TIME = 1

# users that supply a valid API key can post 10 requests per second
ENTREZ_RATE_LOW = 3
ENTREZ_RATE_HIGH = 10

# maximum number of IDs that can be requested in one operation
ENTREZ_MAX_UID = 200


def entrez_request(action, params=None):
    """
    Helper function to ensure that we never exceed the rate limit.
    """
    params = params or dict()

    # tell NCBI which application is making these requests
    params["tool"] = ENTREZ_TOOL
    params["email"] = ENTREZ_EMAIL

    with open(".snakemake/config.yaml") as fin:
        config = yaml.safe_load(fin)

    if config.get("api_key"):
        # append the user specified api_key
        params["api_key"] = config["api_key"]

    url = ENTREZ_URL + action

    if len(params.get("id", [])) > ENTREZ_MAX_UID:
        raise RuntimeError(f"List of Entrez IDs exceeds the maximum: {ENTREZ_MAX_UID}")

    if config.get("debug"):
        # turn into a get request
        get = dict((key, value if isinstance(value, (str, int)) else ",".join(value)) for key, value in params.items())
        print(f"{url}?{urlencode(get)}")

    # make the request
    r = requests.post(url, params)

    if not r.ok:
        r.raise_for_status()

    time.sleep(ENTREZ_WAIT_TIME)

    return r


def entrez_esearch(database, query):
    """
    Execute an Entrez esearch query and return the search keys
    """
    r = entrez_request("esearch.fcgi", {"db": database, "term": query, "usehistory": "y"})

    # parse the XML result
    etree = ElementTree.XML(r.text)

    # get the search keys
    key = etree.find("QueryKey").text
    webenv = etree.find("WebEnv").text
    id_list = [id.text for id in etree.findall(".//Id")]

    return key, webenv, id_list


def entrez_esummary(database, query_key, webenv):
    """
    Fetch the Entrez esummary records for an esearch query.
    """
    r = entrez_request("esummary.fcgi", {"db": database, "query_key": query_key, "WebEnv": webenv})

    return ElementTree.XML(r.text)


def entrez_efetch(database, id_list):
    """
    Fetch the Entrez records for a list of IDs.
    """
    r = entrez_request("efetch.fcgi", {"db": database, "id": id_list})

    return ElementTree.XML(r.text)


def entrez_assembly_ftp(accession, force=False):
    """
    Get an NCBI ftp url from the assembly database.
    """

    # query the assembly database to get the latest assembly for this accession code
    key, webenv, id_list = entrez_esearch("assembly", accession + ' AND "latest"[filter]')

    if len(id_list) > 1:
        # should never happen, but...
        raise RuntimeError(f"Multiple assembly accessions found for '{accession}': {id_list}")

    elif len(id_list) == 0:
        # no entry in the assembly database for this accession code
        return ""

    # fetch the summary record for the assembly
    r = entrez_request("esummary.fcgi", {"db": "assembly", "id": id_list})

    # parse the XML result
    etree = ElementTree.XML(r.text)

    # check if the assembly is anomalous
    anomalous = [reason.text for reason in etree.findall(".//Anomalous/Property")]

    if anomalous:
        message = f"Assembly '{accession}' has been marked as anomalous for reasons: '{'; '.join(anomalous)}'"

        if force:
            print(f"WARNING: {message}", file=sys.stderr)
        else:
            raise RuntimeError(message)

    # preference RefSeq URLs over GenBank URLs
    ftp_stub = etree.find(".//FtpPath_RefSeq") or etree.find(".//FtpPath_GenBank")

    if ftp_stub is None:
        return ""

    # append the fasta filename
    ftp_url = os.path.join(ftp_stub.text, os.path.basename(ftp_stub.text) + "_genomic.fna.gz")

    return ftp_url


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
        raise RuntimeError(f"Could not resolve the accession range '{first}-{last}' for master record '{accession}'")


def entrez_xml_to_dict(etree):
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
