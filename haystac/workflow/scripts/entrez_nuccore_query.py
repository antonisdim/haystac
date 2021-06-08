#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import csv

from haystac.workflow.scripts.entrez_utils import (
    entrez_esearch,
    entrez_esummary,
    entrez_xml_to_dict,
)
from haystac.workflow.scripts.utilities import print_error


def entrez_nuccore_query(query, output_file):
    """
    Query the NCBI nuccore database to get a list of sequence accessions and their metadata.
    """
    # execute the search
    key, webenv, id_list = entrez_esearch("nuccore", query)

    # check that there was at least one record found
    if len(id_list) == 0:
        print_error(f"The --query '{query}' returned no results")

    # fetch the results
    etree = entrez_esummary("nuccore", key, webenv)

    # convert the ElementTree into a a list of dicts
    data = entrez_xml_to_dict(etree)

    with open(output_file, "w") as fout:
        w = csv.DictWriter(fout, data[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(data)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    entrez_nuccore_query(
        query=snakemake.config["query"],
        output_file=snakemake.output[0],
    )
