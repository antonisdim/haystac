#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import time

import pandas as pd

from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import guts_of_entrez


def entrez_download_sequence(accession, config, output_file):
    """
    Fetch the reference genome from NCBI.
    """

    print("The sequences for the database are being selected ...", file=sys.stderr)

    Entrez.email = config['entrez']['email']

    record = guts_of_entrez('nuccore', 'text', [accession], config)

    with open(output_file, 'w') as fout:
        fout.write(record.read())


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_download_sequence(
        accession=snakemake.wildcards.accession,
        config=snakemake.config,
        output_file=snakemake.output[0]
    )
