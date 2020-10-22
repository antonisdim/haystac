#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import gzip
from xml.etree import ElementTree

import os
import sys
import urllib
from Bio import bgzf

from haystack.workflow.scripts.entrez_utils import entrez_esearch, entrez_range_accessions, entrez_request


def entrez_download_sequence(accession, output_file):
    """
    Fetch the Entrez fasta record for a nuccore accession.
    """

    # open the output file stream
    bgzip_fout = bgzf.open(output_file, "wt")

    # query the assembly database to see if there is an FTP url we can use
    key, webenv, id_list = entrez_esearch("assembly", accession + ' AND "latest refseq"[filter]')

    if len(id_list) == 1:
        # get the ID of the corresponding assembly
        assembly_id = id_list.pop()

        # fetch the assembly record
        r = entrez_request(f"esummary.fcgi?db=assembly&id={assembly_id}")

        # parse the XML result
        etree = ElementTree.XML(r.text)

        # preference RefSeq URLs over GenBank URLs
        ftp_stub = etree.find(".//FtpPath_RefSeq").text or etree.find(".//FtpPath_GenBank").text

        if ftp_stub:
            # add missing filename
            ftp_url = os.path.join(ftp_stub, os.path.basename(ftp_stub) + "_genomic.fna.gz")

            # read the FTP stream, unzip the contents and write them one line at a time to our bgzip file
            with gzip.open(urllib.request.urlretrieve(ftp_url)[0]) as fin:
                for line in fin:
                    print(line.strip().decode("utf-8"), file=bgzip_fout)

    else:
        # fetch the fasta record from nuccore
        r = entrez_request(f"efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text")

        # the fasta may be empty if this is a "master record" containing multiple other records (e.g. NZ_APLR00000000.1)
        if len(r.text.strip()) == 0:

            # get the full GenBank XML record
            r = entrez_request(f"efetch.fcgi?db=nuccore&id={accession}&rettype=gb&retmode=xml")

            # parse the XML result
            etree = ElementTree.XML(r.text)

            # get the first and last accession codes for this master record
            first = etree.find(".//GBAltSeqItem_first-accn").text
            last = etree.find(".//GBAltSeqItem_last-accn").text

            if first is None or last is None:
                print(f"ERROR: Could not download the fasta file for {accession}", file=sys.stderr)
                exit(1)

            # get all the related accession codes
            accessions = ",".join(entrez_range_accessions(accession, first, last))

            # fetch all the accessions at once
            r = entrez_request(f"efetch.fcgi?db=nuccore&id={accessions}&rettype=fasta&retmode=text")

        # write the fasta data to our bgzip file
        print(r.text, file=bgzip_fout)

    # close the bgzip file
    bgzip_fout.close()


if __name__ == "__main__":
    entrez_download_sequence(accession=snakemake.wildcards.accession, output_file=snakemake.output[0])
