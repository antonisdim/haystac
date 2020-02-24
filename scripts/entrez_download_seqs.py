#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import time

import pandas as pd

from Bio import Entrez

sys.path.append(os.getcwd())

from scripts.entrez_utils import guts_of_entrez


def make_refseq_database(sequences_file, refseq_database):

    sequences = pd.read_csv(sequences_file, sep='\t')

    for key, seq in sequences.iterrows():
        orgname = seq['TSeq_orgname'].replace(" ", ".")
        accession = seq['TSeq_accver']

        if not os.path.exists(refseq_database + orgname):
            os.makedirs(refseq_database + orgname)

        fetch_refseq_acc(refseq_database, orgname, accession)


def fetch_refseq_acc(path_to_store, folder, acc):
    """Fetch the reference genome from NCBI."""

    config = snakemake.config
    Entrez.email = config['entrez']['email']

    print(acc, file=sys.stderr)
    config['entrez']['retmode'] = 'text'
    records = guts_of_entrez('nuccore', [acc], config)

    name_of_the_file = os.path.join(path_to_store, folder, acc + ".fasta")
    with open(name_of_the_file, 'w') as seqs:
        print(records.read(), file=seqs)
    time.sleep(1)


# todo check with how I can see if the query has been executed before - as Evan wrote as a comment on the google doc
# todo can I individually check every folder that has been dowloaded, in order to sometimes just add missing taxa
#  instead of redowloading them
# todo REFSEQ_DATABASE is something that was previously in the config file. I'll try
#  and move it there in this version too

if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    print("The sequences for the database are being selected ...", file=sys.stderr)

    make_refseq_database(
        sequences_file=snakemake.input[0],
        refseq_database=snakemake.output[0]
    )
