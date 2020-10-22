#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import csv
import os
import pandas as pd
import sys
from Bio import Entrez

from haystack.workflow.scripts.entrez_utils import (
    chunker,
    guts_of_entrez,
    ENTREZ_DB_TAXA,
    ENTREZ_RETMODE_XML,
    ENTREZ_RETTYPE_FASTA,
)


def entrez_taxonomy_query(config, nuccore_file, output_file):
    """
    Query the NCBI taxonomy database to get a taxa details for all nuccore sequences.
    """

    assert os.stat(nuccore_file).st_size, f"The nuccore_query count file is empty {nuccore_file}"

    Entrez.email = config["email"]
    retmax = config["batchsize"]

    # load the unique list of taxa from the nuccore resultset
    accessions = pd.read_csv(nuccore_file, sep="\t", usecols=["TaxId"], squeeze=True).unique()

    with open(output_file, "w") as fout:
        columns = [
            "TSeq_taxid",
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "subspecies",
            "serovar",
        ]
        w = csv.DictWriter(fout, columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        resultset = len(accessions)

        for chunk in chunker(accessions, retmax):
            print(
                f"Remaining sequences to have their taxids downloaded {resultset}\n", file=sys.stderr,
            )
            records = guts_of_entrez(ENTREZ_DB_TAXA, ENTREZ_RETMODE_XML, ENTREZ_RETTYPE_FASTA, chunk, retmax)

            for node in records:
                taxon = dict()
                taxon["TSeq_taxid"] = node["TaxId"]
                for item in node["LineageEx"]:
                    if item["Rank"] in columns:
                        taxon[item["Rank"]] = item["ScientificName"]
                    elif item["Rank"] == "no rank":
                        taxon["serovar"] = item["ScientificName"]

                if node["Rank"] in columns:
                    taxon[node["Rank"]] = node["ScientificName"]

                if "serovar" not in taxon.keys():
                    if node["Rank"] == "no rank":
                        taxon["serovar"] = node["ScientificName"]

                if taxon["species"] and not taxon.get("subspecies"):
                    taxon["subspecies"] = taxon["species"] + " ssp."

                if taxon["species"] and not taxon.get("serovar"):
                    taxon["serovar"] = ""

                w.writerow(taxon)

            resultset -= retmax
            print("done for this slice\n", file=sys.stderr)

    print("COMPLETE", file=sys.stderr)


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    entrez_taxonomy_query(
        config=snakemake.config, nuccore_file=snakemake.input[0], output_file=snakemake.output[0],
    )
