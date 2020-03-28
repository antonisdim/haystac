#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import urllib.error
import time

from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import guts_of_entrez, ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_TEXT, ENTREZ_RETTYPE_FASTA

TOO_MANY_REQUESTS_WAIT = 5
MAX_RETRY_ATTEMPTS = 3


def entrez_download_sequence(accession, config, output_file, attempt=1):
    """
    Fetch the reference genome from NCBI.
    """
    print("The sequences for the database are being selected ...", file=sys.stderr)

    Entrez.email = config['entrez']['email']

    try:
        records = guts_of_entrez(ENTREZ_DB_NUCCORE, ENTREZ_RETMODE_TEXT, ENTREZ_RETTYPE_FASTA, [accession],
                                 batch_size=1)
        with gzip.open(output_file, 'wt') as fout:
            for fasta in records:
                fout.write(fasta)

    except(urllib.error.HTTPError, e):
        if e.code == 429:

            attempt += 1

            if attempt > MAX_RETRY_ATTEMPTS:
                print("Exceeded maximum attempts {}...".format(attempt), file=sys.stderr)
                return None
            else:
                time.sleep(TOO_MANY_REQUESTS_WAIT)
                entrez_download_sequence(accession, config, output_file, attempt)

        else:
            raise RuntimeError("There was a urllib.error.HTTPError with code {}".format(e))



if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_download_sequence(
        accession=snakemake.wildcards.accession,
        config=snakemake.config,
        output_file=snakemake.output[0]
    )
