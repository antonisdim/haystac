#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from haystac.workflow.scripts.utilities import REGEX_BLACKLIST


def entrez_pick_sequences(config, nuccore_file, taxa_file, output_file):
    accessions = pd.read_csv(nuccore_file, sep="\t")
    taxa = pd.read_csv(taxa_file, sep="\t")
    rank = config["rank"]

    sequences = pd.merge(accessions, taxa, on="TaxId", how="outer")

    sequences = sequences[~sequences[rank].isnull()]

    selected_sequences = sequences.loc[
        sequences.groupby(rank)["Length"].idxmax(),
        ["species", "AccessionVersion"],
    ]

    if config["refseq_rep"]:
        refseq_genomes = pd.read_csv(config["db_output"] + "/entrez/refseq-genomes.tsv", sep="\t")
        genbank_genomes = pd.read_csv(config["db_output"] + "/entrez/genbank-genomes.tsv", sep="\t")
        assemblies = pd.read_csv(config["db_output"] + "/entrez/assemblies.tsv", sep="\t")
        refseq_plasmids = pd.read_csv(config["db_output"] + "/entrez/refseq-plasmids.tsv", sep="\t")
        genbank_plasmids = pd.read_csv(config["db_output"] + "/entrez/genbank-plasmids.tsv", sep="\t")

        # the entrez query might give a different accession for a certain species than the refseq rep one and
        # I don't want that. If the species exists in the refseq I want to keep the refseq records
        selected_sequences = selected_sequences[
            ~selected_sequences["species"].replace(REGEX_BLACKLIST, "_", regex=True).isin(refseq_genomes.species)
        ]
        selected_sequences = selected_sequences[
            ~selected_sequences["species"].replace(REGEX_BLACKLIST, "_", regex=True).isin(genbank_genomes.species)
        ]
        selected_sequences = selected_sequences[
            ~selected_sequences["species"].replace(REGEX_BLACKLIST, "_", regex=True).isin(assemblies.species)
        ]
        selected_sequences = selected_sequences[
            ~selected_sequences["species"].replace(REGEX_BLACKLIST, "_", regex=True).isin(refseq_plasmids.species)
        ]
        selected_sequences = selected_sequences[
            ~selected_sequences["species"].replace(REGEX_BLACKLIST, "_", regex=True).isin(genbank_plasmids.species)
        ]

    selected_sequences["species"] = selected_sequences["species"].replace(REGEX_BLACKLIST, "_", regex=True)

    if config["sequences"]:
        custom_fasta_paths = pd.read_csv(
            config["sequences"],
            sep="\t",
            header=None,
            names=["species", "accession", "path"],
        )

        selected_sequences = selected_sequences[(~selected_sequences["species"].isin(custom_fasta_paths["species"]))]

    if config["accessions"]:
        custom_accessions = pd.read_csv(
            config["accessions"],
            sep="\t",
            header=None,
            names=["species", "accession"],
        )

        selected_sequences = selected_sequences[(~selected_sequences["species"].isin(custom_accessions["species"]))]

    selected_sequences.to_csv(output_file, sep="\t", header=True, index=False)


if __name__ == "__main__":
    # noinspection PyUnresolvedReferences
    entrez_pick_sequences(
        config=snakemake.config,
        nuccore_file=snakemake.input[0],
        taxa_file=snakemake.input[1],
        output_file=snakemake.output[0],
    )
