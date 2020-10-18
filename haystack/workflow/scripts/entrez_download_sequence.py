#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import argparse
from xml.etree import ElementTree

import requests

# base url of the Entrez web service
ENTREZ_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def entrez_range_stings(first, last):
    """
    Return a range between two strings (e.g. ABC001 -> ABC005)
    """
    # sanity check that the items are equal length and in order
    assert len(first) == len(last) and first < last

    # find the index of the first difference
    idx = [i for i in range(len(first)) if first[i] != last[i]][0]

    # return the range
    return [f"{first[:idx]}{item}" for item in range(int(first[idx:]), int(last[idx:]) + 1)]


def entrez_download_sequence(database, accession):
    """
    Fetch the Entrez fasta record for a nuccore accession.
    """
    r = requests.get(ENTREZ_URL + f"efetch.fcgi?db={database}&id={accession}&rettype=fasta&retmode=text")

    fasta = r.text

    # result can be empty if this is a "master record" containing multiple other records (e.g. NZ_APLR00000000.1)
    if len(fasta.strip()) == 0:
        # get the full GenBank XML record
        r = requests.get(ENTREZ_URL + f"efetch.fcgi?db={database}&id={accession}&rettype=gb&retmode=xml")

        # parse the XML result
        etree = ElementTree.XML(r.text)

        # get the first and last accession codes for this master record
        first = etree.find(".//GBAltSeqItem_first-accn").text
        last = etree.find(".//GBAltSeqItem_last-accn").text

        for accession in entrez_range_stings(first, last):
            r = requests.get(ENTREZ_URL + f"efetch.fcgi?db={database}&id={accession}&rettype=fasta&retmode=text")
            print(r.text)
    else:
        print(r.text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download seq data from Entrez.")
    parser.add_argument("-d", "--database", help="NCBI database")
    parser.add_argument("-a", "--accession", help="Accession to be downloaded")
    args = parser.parse_args()

    entrez_download_sequence(args.database, args.accession)
