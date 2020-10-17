#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from Bio import Entrez

from haystack.workflow.scripts.entrez_utils import (
    ENTREZ_DB_NUCCORE,
    ENTREZ_RETMODE_XML,
    ENTREZ_RETTYPE_GB,
    ENTREZ_RETMAX,
)


def entrez_find_accessions(config, output_file):
    """Function to find all the accessions related to our NCBI query"""

    Entrez.email = config["email"]

    entrez_query = config["query"]

    handle = Entrez.esearch(
        db=ENTREZ_DB_NUCCORE,
        term=entrez_query,
        retmax=ENTREZ_RETMAX,
        idtype="acc",
        rettype=ENTREZ_RETTYPE_GB,
        retmode=ENTREZ_RETMODE_XML,
    )

    handle_reader = Entrez.read(handle)

    accessions = handle_reader["IdList"]

    with open(output_file, "w") as fout:
        acc_col = ["GBSeq_accession-version"]
        acc_col.extend(accessions)

        print("\n".join(acc_col), file=fout)


if __name__ == "__main__":
    entrez_find_accessions(config=snakemake.config, output_file=snakemake.output[0])
