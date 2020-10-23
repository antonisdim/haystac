#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import csv
import pandas as pd

from haystack.workflow.scripts.entrez_utils import entrez_esearch


def entrez_invalid_assemblies(assemblies, output):
    assemblies_file = pd.read_csv(assemblies, sep="\t")

    with open(output, "w") as fout:
        columns = ["species", "AccessionVersion"]
        w = csv.DictWriter(fout, columns, delimiter="\t")
        w.writeheader()

        # TODO why do this in a loop? much quicker to post all the IDs at ounce
        for key, acc in assemblies_file.iterrows():
            # query the assembly database to confirm that this accession is still valid
            _, _, id_list = entrez_esearch("assembly", acc["AccessionVersion"] + ' AND "latest refseq"[filter]')

            if len(id_list) == 0:
                row = dict()
                row["species"] = acc["species"]
                row["AccessionVersion"] = acc["AccessionVersion"]
                w.writerow(row)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    entrez_invalid_assemblies(assemblies=snakemake.input[0], output=snakemake.output[0])
