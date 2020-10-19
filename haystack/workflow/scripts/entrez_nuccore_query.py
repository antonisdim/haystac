#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import csv
import sys

from haystack.workflow.scripts.entrez_utils import entrez_esearch, entrez_esummary


def element_tree_to_dict(etree):
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


def entrez_nuccore_query(query, output_file):
    """
    Query the NCBI nuccore database to get a list of sequence accessions and their metadata.
    """
    # execute the search
    key, webenv, _ = entrez_esearch("nuccore", query)

    # fetch the results
    etree = entrez_esummary("nuccore", key, webenv)

    # convert the ElementTree into a a list of dicts
    data = element_tree_to_dict(etree)

    with open(output_file, "w") as fout:
        w = csv.DictWriter(fout, data[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(data)


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    entrez_nuccore_query(
        query=snakemake.config["query"], output_file=snakemake.output[0],
    )
