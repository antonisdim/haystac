#!/usr/bin/env python
# -*- coding: utf-8 -*

import csv
import os
import sys

import pandas as pd
from Bio import Entrez

sys.path.append(os.getcwd())


from scripts.entrez_utils import ENTREZ_DB_ASSEMBLY

TOO_MANY_REQUESTS_WAIT = 5
MAX_RETRY_ATTEMPTS = 3


def entrez_invalid_assemblies(config, assemblies, output):
    assemblies_file = pd.read_csv(assemblies, sep='\t')

    Entrez.email = config['entrez_email']

    with open(output, 'w') as fout:
        columns = ['species', 'GBSeq_accession-version']
        w = csv.DictWriter(fout, columns, delimiter='\t', extrasaction="ignore")
        w.writeheader()

        for key, acc in assemblies_file.iterrows():

            invalid_assemblies = dict()

            handle = Entrez.esearch(db=ENTREZ_DB_ASSEMBLY,
                term=acc['GBSeq_accession-version'] + ' AND "latest refseq"[filter]')
            assembly_record = Entrez.read(handle)

            if not len(assembly_record['IdList']) > 0:
                invalid_assemblies['species'] = acc['species']
                invalid_assemblies['GBSeq_accession-version'] = acc['GBSeq_accession-version']

                print(invalid_assemblies, file=sys.stderr)
                w.writerow(invalid_assemblies)


if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    entrez_invalid_assemblies(
        config=snakemake.config,
        assemblies=snakemake.input[0],
        output=snakemake.output[0]
    )
