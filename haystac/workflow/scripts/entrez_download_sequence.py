#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import gzip
import urllib.error
import urllib.request
from xml.etree import ElementTree

import requests
from Bio import bgzf

from haystac.workflow.scripts.entrez_utils import (
    entrez_range_accessions,
    entrez_request,
    entrez_assembly_ftp,
    entrez_find_replacement_accession,
    ENTREZ_MAX_UID,
    ENTREZ_MAX_ATTEMPTS,
)
from haystac.workflow.scripts.utilities import chunker, print_error


def download_entrez_ftp(ftp_url, output_file, attempt=1):
    """
    Read the FTP stream, unzip the contents and write them one line at a time to our bgzip file
    """
    try:
        # download the assembly file to a temp location
        temp_file, _ = urllib.request.urlretrieve(ftp_url)

        # stream the gzip file and recode as bgzip
        with gzip.open(temp_file) as fin, bgzf.open(output_file, "w") as fout:
            for line in fin:
                print(line.strip().decode("utf-8"), file=fout)

    except urllib.error.URLError as error:
        if attempt < ENTREZ_MAX_ATTEMPTS:
            # try downloading it again
            download_entrez_ftp(ftp_url, output_file, attempt + 1)
        else:
            raise error

    finally:
        # delete any temp files which are left behind by urlretrieve, and may fill up the `/tmp` folder
        urllib.request.urlcleanup()


def entrez_download_sequence(accession, output_file, force=False, mtdna=False):
    """
    Fetch the Entrez fasta record for a nuccore accession.
    """

    # query the assembly database to see if there is an FTP url we can use
    ftp_url = entrez_assembly_ftp(accession, force) if not mtdna else ""

    try:
        if ftp_url:
            download_entrez_ftp(ftp_url, output_file)
            return

    except urllib.error.URLError:
        pass

    try:
        # fetch the fasta record from nuccore
        r = entrez_request("efetch.fcgi", {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"})

        # the fasta may be empty if this is a "master record" containing multiple other records (NZ_APLR00000000.1)
        if len(r.text.strip()) > 1:
            with bgzf.open(output_file, "w") as fout:
                print(r.text, file=fout)

            return

    except requests.exceptions.HTTPError:
        pass

    try:
        # get the full GenBank XML record
        r = entrez_request("efetch.fcgi", {"db": "nuccore", "id": accession, "rettype": "gb", "retmode": "xml"})

    except requests.exceptions.HTTPError:
        # check for a replacement accession (there may be a newer version if this a WGS project)
        updated_accession = entrez_find_replacement_accession(accession)

        # download the updated accession instead
        entrez_download_sequence(updated_accession, output_file, force)

        return

    # parse the XML result
    etree = ElementTree.XML(r.text)

    # get the first and last accession codes for this master record
    first = etree.find(".//GBAltSeqItem_first-accn")
    last = etree.find(".//GBAltSeqItem_last-accn")

    if first is None or last is None:
        print_error(f"Could not download the fasta file for {accession}")

    # get all the related accession codes
    accessions = entrez_range_accessions(accession, first.text, last.text)

    try:
        with bgzf.open(output_file, "w") as fout:
            # fetch all the accessions in batches
            for id_list in chunker(accessions, ENTREZ_MAX_UID):
                r = entrez_request(
                    "efetch.fcgi", {"db": "nuccore", "id": id_list, "rettype": "fasta", "retmode": "text"},
                )

                # write the fasta data to our bgzip file
                print(r.text, file=fout)

    except requests.exceptions.HTTPError:
        print_error(
            f"Could not download the accession range '{first.text}-{last.text}' for master record '{accession}'"
        )


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    entrez_download_sequence(
        accession=snakemake.wildcards.accession,
        output_file=snakemake.output[0],
        force=snakemake.config["force_accessions"],
        mtdna=snakemake.config["mtDNA"],
    )
