#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import gzip
from xml.etree import ElementTree

import requests
import sys
from Bio import bgzf
from urllib.request import urlretrieve

from haystack.workflow.scripts.entrez_utils import (
    entrez_range_accessions,
    entrez_request,
    entrez_assembly_ftp,
    ENTREZ_MAX_UID,
)
from haystack.workflow.scripts.utilities import chunker


def entrez_download_sequence(accession, output_file, force=False):
    """
    Fetch the Entrez fasta record for a nuccore accession.
    """

    # open the output file stream
    with bgzf.open(output_file, "wt") as bgzip_fout:

        # query the assembly database to see if there is an FTP url we can use
        ftp_url = entrez_assembly_ftp(accession, force)

        if ftp_url:
            # read the FTP stream, unzip the contents and write them one line at a time to our bgzip file
            with gzip.open(urlretrieve(ftp_url)[0]) as fin:
                for line in fin:
                    print(line.strip().decode("utf-8"), file=bgzip_fout)

        else:
            try:
                # fetch the fasta record from nuccore
                r = entrez_request(
                    "efetch.fcgi", {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"}
                )
                fasta = r.text
            except requests.exceptions.HTTPError:
                fasta = ""

            # the fasta may be empty if this is a "master record" containing multiple other records (NZ_APLR00000000.1)
            if len(fasta.strip()) == 0:

                # get the full GenBank XML record
                r = entrez_request("efetch.fcgi", {"db": "nuccore", "id": accession, "rettype": "gb", "retmode": "xml"})

                # parse the XML result
                etree = ElementTree.XML(r.text)

                # get the first and last accession codes for this master record
                first = etree.find(".//GBAltSeqItem_first-accn")
                last = etree.find(".//GBAltSeqItem_last-accn")

                if first is None or last is None:
                    raise RuntimeError(f"ERROR: Could not download the fasta file for {accession}")

                # get all the related accession codes
                accessions = entrez_range_accessions(accession, first.text, last.text)

                try:
                    # fetch all the accessions in batches
                    for id_list in chunker(accessions, ENTREZ_MAX_UID):
                        r = entrez_request(
                            "efetch.fcgi", {"db": "nuccore", "id": id_list, "rettype": "fasta", "retmode": "text"}
                        )
                        fasta += r.text

                except requests.exceptions.HTTPError:
                    raise RuntimeError(
                        f"Could not download the accession range '{first.text}-{last.text}' "
                        f"for master record '{accession}'"
                    )

            # write the fasta data to our bgzip file
            print(fasta, file=bgzip_fout)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    entrez_download_sequence(
        accession=snakemake.wildcards.accession,
        output_file=snakemake.output[0],
        force=snakemake.config["force_accessions"],
    )
