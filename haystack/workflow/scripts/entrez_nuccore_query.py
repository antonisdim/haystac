#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import csv

from haystack.workflow.scripts.entrez_utils import entrez_esearch, entrez_esummary_webenv, entrez_xml_to_dict


def entrez_nuccore_query(query, output_file):
    """
    Query the NCBI nuccore database to get a list of sequence accessions and their metadata.
    """
    # execute the search
    key, webenv, _ = entrez_esearch("nuccore", query)

    # fetch the results
    etree = entrez_esummary_webenv("nuccore", key, webenv)

    # convert the ElementTree into a a list of dicts
    data = entrez_xml_to_dict(etree)

    with open(output_file, "w") as fout:
        w = csv.DictWriter(fout, data[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(data)


if __name__ == "__main__":
    entrez_nuccore_query(
        query=snakemake.config["query"], output_file=snakemake.output[0],
    )
