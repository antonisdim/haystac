#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import guts_of_entrez, ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_TEXT, ENTREZ_RETTYPE_FASTA


def entrez_download_sequence(accession, config, output_file):
    """
    Fetch the reference genome from NCBI.
    """
    print("The sequences for the database are being selected ...", file=sys.stderr)

    Entrez.email = config['entrez']['email']

    records = guts_of_entrez(ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_TEXT, ENTREZ_RETTYPE_FASTA, [accession], batch_size=1)

    with open(output_file, 'w') as fout:
        for fasta in records:
            fout.write(fasta)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_download_sequence(
        accession=snakemake.wildcards.accession,
        config=snakemake.config,
        output_file=snakemake.output[0]
    )
