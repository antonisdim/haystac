#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import csv
import os

import pandas as pd

from haystac.workflow.scripts.entrez_utils import entrez_efetch, ENTREZ_MAX_UID
from haystac.workflow.scripts.utilities import chunker


def entrez_taxonomy_query(nuccore_file, output_file):
    """
    Query the NCBI taxonomy database to get the taxa details for all the nuccore sequences.
    """

    assert os.stat(nuccore_file).st_size, f"The nuccore_query count file is empty {nuccore_file}"

    # load the unique list of taxa from the nuccore resultset
    accessions = [tax_id for tax_id in pd.read_csv(nuccore_file, sep="\t", usecols=["TaxId"], squeeze=True).unique()]

    with open(output_file, "w") as fout:
        columns = [
            "TaxId",
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

        for id_list in chunker(accessions, ENTREZ_MAX_UID):

            # fetch the taxonomy metadata
            etree = entrez_efetch("taxonomy", id_list)

            for taxon in etree:
                row = dict()
                row["TaxId"] = taxon.find("TaxId").text

                # extract the taxonomic hierarchy
                for lineage in taxon.find("LineageEx"):
                    lineage_rank = lineage.find("Rank").text
                    if lineage_rank in columns:
                        row[lineage_rank] = lineage.find("ScientificName").text
                    elif lineage_rank == "no rank":
                        row["serovar"] = lineage.find("ScientificName").text

                # get the taxon rank and name
                taxon_rank = taxon.find("Rank").text
                taxon_name = taxon.find("ScientificName").text

                if taxon_rank in columns:
                    row[taxon_rank] = taxon_name

                if "serovar" not in row.keys():
                    if taxon_rank == "no rank":
                        row["serovar"] = taxon_name

                if row["species"] and not row.get("subspecies"):
                    row["subspecies"] = row["species"] + " ssp."

                if row["species"] and not row.get("serovar"):
                    row["serovar"] = ""

                w.writerow(row)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    entrez_taxonomy_query(
        nuccore_file=snakemake.input[0],
        output_file=snakemake.output[0],
    )
