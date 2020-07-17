#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import csv
import os
import sys

from Bio import Entrez

sys.path.append(os.getcwd())
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from scripts.entrez_utils import (
    guts_of_entrez,
    ENTREZ_DB_NUCCORE,
    ENTREZ_RETMODE_XML,
    ENTREZ_RETTYPE_GB,
    ENTREZ_RETMAX,
)


def entrez_find_accessions(config, output_file):
    Entrez.email = config["entrez_email"]

    entrez_query = config["entrez_query"]

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
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    entrez_find_accessions(config=snakemake.config, output_file=snakemake.output[0])
