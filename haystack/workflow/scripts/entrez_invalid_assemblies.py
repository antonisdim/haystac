#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import csv
import os
import sys
import time
import urllib.error

import pandas as pd
from Bio import Entrez

sys.path.append(os.getcwd())
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from scripts.entrez_utils import ENTREZ_DB_ASSEMBLY

TOO_MANY_REQUESTS_WAIT = 10
MAX_RETRY_ATTEMPTS = 5


def entrez_invalid_assemblies(config, assemblies, output, attempt=1):
    assemblies_file = pd.read_csv(assemblies, sep="\t")

    Entrez.email = config["email"]

    with open(output, "w") as fout:
        columns = ["species", "GBSeq_accession-version"]
        w = csv.DictWriter(fout, columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        for key, acc in assemblies_file.iterrows():

            invalid_assemblies = dict()

            try:
                handle = Entrez.esearch(
                    db=ENTREZ_DB_ASSEMBLY, term=acc["GBSeq_accession-version"] + ' AND "latest refseq"[filter]',
                )
                assembly_record = Entrez.read(handle)

                if not len(assembly_record["IdList"]) > 0:
                    invalid_assemblies["species"] = acc["species"]
                    invalid_assemblies["GBSeq_accession-version"] = acc["GBSeq_accession-version"]

                    print(invalid_assemblies, file=sys.stderr)
                    w.writerow(invalid_assemblies)

            except urllib.error.HTTPError as e:
                if e.code == 429:

                    attempt += 1

                    if attempt > MAX_RETRY_ATTEMPTS:
                        print(
                            "Exceeded maximum attempts {}...".format(attempt), file=sys.stderr,
                        )
                        return None
                    else:
                        time.sleep(TOO_MANY_REQUESTS_WAIT)
                        entrez_invalid_assemblies(config, assemblies, output, attempt)

                else:
                    raise RuntimeError("There was a urllib.error.HTTPError with code {}".format(e))


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    entrez_invalid_assemblies(
        config=snakemake.config, assemblies=snakemake.input[0], output=snakemake.output[0],
    )