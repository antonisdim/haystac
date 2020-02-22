#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pandas as pd
import sys
import time
from Bio import Entrez

sys.path.append(os.getcwd())
from scripts.entrez_utils import chunker, guts_of_entrez

def make_refseq_database(input_file):
    # make the folder that will contain the taxonomic info
    # if not os.path.exists(REFSEQ_DATABASE):
    #     os.makedirs(REFSEQ_DATABASE)

    refseq_database = snakemake.output[0]

    result = {}
    with open(input_file) as selected:
        for line in selected:
            line = line.rstrip()
            (key, val) = line.split('\t')
            result[key] = val

    for species, accessions in result.items():
        species = species.replace(" ", ".")

        if not os.path.exists(refseq_database + species):
            os.makedirs(refseq_database + species)

        fetch_refseq_acc(refseq_database, species, accessions)


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
#  todo can I individually check every folder that has been dowloaded, in order to sometimes just add missing taxa
#   instead of redowloading them todo REFSEQ_DATABASE is something that was previously in the config file. I'll try
#    and move it there in this version too

if __name__ == '__main__':
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], 'w')

    print("The sequences for the database are being selected ...", file=sys.stderr)

    make_refseq_database(
        input_file=snakemake.input[0]
    )
