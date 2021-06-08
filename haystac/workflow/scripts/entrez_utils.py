#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
import time
from urllib.parse import urlencode
from xml.etree import ElementTree

import requests

from haystac.workflow.scripts.utilities import print_warning, print_error, get_smk_config

# base url of the Entrez web service
ENTREZ_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# tool registration with NCBI, see https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Frequency_Timing_and_Registrati
ENTREZ_TOOL = "haystac"
ENTREZ_EMAIL = "antonisdim41@gmail.com"

# maximum 1 request per X seconds
ENTREZ_WAIT_TIME = 1.1

# users that supply a valid API key can post 10 requests per second
ENTREZ_RATE_LOW = 3
ENTREZ_RATE_HIGH = 10

# maximum number of IDs that can be requested in one operation
ENTREZ_MAX_UID = 200

# maximum number of retries when we encounter NCBI server/network error
ENTREZ_MAX_ATTEMPTS = 2

# check if the error is because of NCBI [101=NewConnectionError, 104=ConnectionResetError, 110=TimeoutError,
# 111=ConnectionRefusedError, 500=InternalServerError]
ENTREZ_ERRORS = [101, 104, 110, 111, 500]


def entrez_request(action, params=None, attempt=1):
    """
    Helper function to ensure that we never exceed the rate limit.
    """
    params = params or dict()

    # tell NCBI which application is making these requests
    params["tool"] = ENTREZ_TOOL
    params["email"] = ENTREZ_EMAIL

    config = get_smk_config()

    if config.get("api_key"):
        # append the user specified api_key
        params["api_key"] = config["api_key"]

    url = ENTREZ_URL + action

    if len(params.get("id", [])) > ENTREZ_MAX_UID:
        print_error(f"List of Entrez IDs exceeds the maximum: {ENTREZ_MAX_UID}")

    if config.get("debug"):
        # turn into a get request
        print(params.items())
        get = dict(
            (key, value if isinstance(value, (str, int)) else ",".join(str(val) for val in value))
            for key, value in params.items()
        )
        print(f"{url}?{urlencode(get)}")

    # make the request
    r = requests.post(url, params)

    # enforce the rate limit (even when the request failed)
    time.sleep(ENTREZ_WAIT_TIME)

    if not r.ok:
        if r.status_code in ENTREZ_ERRORS and attempt < ENTREZ_MAX_ATTEMPTS:
            return entrez_request(action, params, attempt + 1)
        else:
            r.raise_for_status()

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
    id_list = [id_node.text for id_node in etree.findall(".//Id")]

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

    # find out if we are looking for a virus to apply the correct filter
    filter_condition = ' AND "latest"[filter] NOT suppressed*'

    config = get_smk_config()

    if config.get("refseq_rep") and config["refseq_rep"] == "viruses":
        # append the virus filter
        filter_condition += " AND viruses[filter] "

    # query the assembly database to get the latest assembly for this accession code
    key, webenv, id_list = entrez_esearch(
        "assembly",
        accession + filter_condition,
    )

    if len(id_list) > 1:
        # should never happen, but...
        msg = f"Multiple assembly accessions found for '{accession}': {id_list}. "

        # if force-accessions is true pick the largest int value, assuming it is also the altest
        if force:
            msg += f"Using assembly: id pair '{accession}': '{max([int(id_num) for id_num in id_list])}'"
            id_list = [str(max([int(id_num) for id_num in id_list]))]
            print_warning(msg)

        # if not raise an error
        else:
            msg += (
                f"Either consider using the `--force-accessions` flag for the largest ID to be picked, "
                f"or the `--exclude-accessions` flag to remove accession '{accession}' from this query."
            )
            print_error(msg)

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
        message = (
            f"Assembly '{accession}' has been marked as anomalous for the following reasons: "
            f"'{'; '.join(anomalous)}'"
        )

        if force:
            print_warning(message)
        else:
            print_error(message)

    refseq = etree.find(".//FtpPath_RefSeq")
    genbank = etree.find(".//FtpPath_GenBank")

    # preference RefSeq URLs over GenBank URLs
    if refseq is not None and refseq.text not in ["", None]:
        ftp_stub = refseq.text
    elif genbank is not None and genbank.text not in ["", None]:
        ftp_stub = genbank.text
    else:
        return ""

    # append the fasta filename
    ftp_url = os.path.join(ftp_stub, os.path.basename(ftp_stub) + "_genomic.fna.gz")

    return ftp_url


def entrez_range_accessions(accession, first, last):
    """
    Return a range between two accession codes (e.g. ABC001 -> ABC005)
    """
    # sanity check that the items are equal length and in order
    assert len(first) == len(last) and first < last

    # find the index of the first difference
    idx = [i for i in range(len(first)) if first[i] != last[i]][0]
    pad = len(first) - idx

    try:
        # return the range
        return [f"{first[:idx]}{str(item).zfill(pad)}" for item in range(int(first[idx:]), int(last[idx:]) + 1)]
    except ValueError:
        print_error(f"Could not resolve the accession range '{first}-{last}' for master record '{accession}'")


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


def entrez_find_replacement_accession(accession):
    """
    If the updated version of an assembly accession if possible
    """

    # try getting a new accession for the master record
    accession_new = accession[:-1] + "0000000"

    try:
        # send a request and see if we get back an xml result
        r = entrez_request(
            "efetch.fcgi",
            {"db": "nuccore", "id": accession_new, "rettype": "gb", "retmode": "xml"},
        )

    except requests.exceptions.HTTPError:
        print_error(f"Could not find either the GenBank record for '{accession}' or an alternative accession")

    # noinspection PyUnboundLocalVariable
    etree = ElementTree.XML(r.text)
    replacement = etree.find(".//GBSeq_accession-version")

    # check that the new accession is a WGS project
    keywords = [keyword.text.lower() for keyword in etree.findall(".//GBKeyword")]

    if replacement is not None and "wgs" in keywords:
        print_warning(f"Replacing the superseded WGS accession '{accession}' with '{replacement.text}'")
        return replacement.text
    else:
        print_error(
            f"Could not find either the GenBank record for '{accession}' or a valid alternative accession. "
            f"Please consider using the `--exclude-accessions` flag to remove accession '{accession}' from this query."
        )
