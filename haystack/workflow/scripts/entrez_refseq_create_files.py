#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import sys

import pandas as pd
import os

sys.path.append(os.getcwd())
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from scripts.rip_utilities import get_accession_ftp_path


def entrez_refseq_create_files(
    config,
    input_file,
    nuccore_genomes_out,
    genbank_genomes_out,
    assemblies_out,
    nuccore_plasmids_out,
    genbank_plasmids_out,
):
    prok_refseq_rep = pd.read_csv(input_file, sep="\t")

    prok_refseq_rep_rmdup = prok_refseq_rep[~prok_refseq_rep["#Species/genus"].duplicated()]

    assemblies = prok_refseq_rep_rmdup.loc[prok_refseq_rep_rmdup["WGS"].notna(), ["#Species/genus", "WGS"]]

    assemblies["#Species/genus"] = assemblies["#Species/genus"].str.replace(" ", "_")
    assemblies["#Species/genus"] = assemblies["#Species/genus"].str.replace("/", "_")
    assemblies["#Species/genus"] = assemblies["#Species/genus"].str.replace("'", "")
    assemblies["#Species/genus"] = assemblies["#Species/genus"].str.replace("(", "")
    assemblies["#Species/genus"] = assemblies["#Species/genus"].str.replace(")", "")

    nuccore = prok_refseq_rep_rmdup.loc[
        prok_refseq_rep_rmdup["Chromosome RefSeq"].notna(), ["#Species/genus", "Chromosome RefSeq"],
    ]

    nuccore["#Species/genus"] = nuccore["#Species/genus"].str.replace(" ", "_")
    nuccore["#Species/genus"] = nuccore["#Species/genus"].str.replace("'", "")
    nuccore["#Species/genus"] = nuccore["#Species/genus"].str.replace("(", "")
    nuccore["#Species/genus"] = nuccore["#Species/genus"].str.replace(")", "")

    genbank = prok_refseq_rep_rmdup.loc[
        prok_refseq_rep_rmdup["Chromosome GenBank"].notna(), ["#Species/genus", "Chromosome GenBank"],
    ]

    genbank["#Species/genus"] = genbank["#Species/genus"].str.replace(" ", "_")
    genbank["#Species/genus"] = genbank["#Species/genus"].str.replace("'", "")
    genbank["#Species/genus"] = genbank["#Species/genus"].str.replace("(", "")
    genbank["#Species/genus"] = genbank["#Species/genus"].str.replace(")", "")

    genbank_filtered = genbank[
        (~genbank["#Species/genus"].isin(assemblies["#Species/genus"]))
        & (~genbank["#Species/genus"].isin(nuccore["#Species/genus"]))
    ].copy()
    genbank_filtered.loc[:, "Chromosome GenBank"] = genbank_filtered["Chromosome GenBank"].str.split(",")
    genbank_exploded = genbank_filtered.explode("Chromosome GenBank")
    genbank_exploded = genbank_exploded[~genbank_exploded["#Species/genus"].duplicated()]

    nuccore_filtered = nuccore[
        (~nuccore["#Species/genus"].isin(assemblies["#Species/genus"]))
        & (~nuccore["#Species/genus"].isin(genbank_filtered["#Species/genus"]))
    ].copy()
    nuccore_filtered.loc[:, "Chromosome RefSeq"] = nuccore_filtered["Chromosome RefSeq"].str.split(",")
    nuccore_exploded = nuccore_filtered.explode("Chromosome RefSeq")
    nuccore_exploded = nuccore_exploded[~nuccore_exploded["#Species/genus"].duplicated()]

    assemblies_filtered = assemblies[
        (~assemblies["#Species/genus"].isin(nuccore_filtered["#Species/genus"]))
        & (~assemblies["#Species/genus"].isin(genbank_filtered["#Species/genus"]))
    ]
    assemblies_filtered.loc[:, "WGS"] = assemblies_filtered["WGS"].str.split(",")
    assemblies_exploded = assemblies_filtered.explode("WGS")
    assemblies_exploded = assemblies_exploded[~assemblies_exploded["#Species/genus"].duplicated()]

    nuccore_plasmids = prok_refseq_rep_rmdup[
        prok_refseq_rep_rmdup["Plasmid RefSeq"].notna() & prok_refseq_rep_rmdup["WGS"].isna()
    ].loc[:, ["#Species/genus", "Plasmid RefSeq"]]

    genbank_plasmids = prok_refseq_rep_rmdup[
        prok_refseq_rep_rmdup["Plasmid GenBank"].notna() & prok_refseq_rep_rmdup["WGS"].isna()
    ].loc[:, ["#Species/genus", "Plasmid GenBank"]]
    genbank_plasmids_filtered = genbank_plasmids[
        ~genbank_plasmids["#Species/genus"].isin(nuccore_plasmids["#Species/genus"])
    ].copy()

    nuccore_plasmids.loc[:, "Plasmid RefSeq"] = nuccore_plasmids["Plasmid RefSeq"].str.split(",")
    nuccore_plasmids_exploded = nuccore_plasmids.explode("Plasmid RefSeq")

    nuccore_plasmids_exploded["#Species/genus"] = nuccore_plasmids_exploded["#Species/genus"].str.replace(" ", "_")
    nuccore_plasmids_exploded["#Species/genus"] = nuccore_plasmids_exploded["#Species/genus"].str.replace("/", "_")
    nuccore_plasmids_exploded["#Species/genus"] = nuccore_plasmids_exploded["#Species/genus"].str.replace("'", "")
    nuccore_plasmids_exploded["#Species/genus"] = nuccore_plasmids_exploded["#Species/genus"].str.replace("(", "")
    nuccore_plasmids_exploded["#Species/genus"] = nuccore_plasmids_exploded["#Species/genus"].str.replace(")", "")

    genbank_plasmids_filtered.loc[:, "Plasmid GenBank"] = genbank_plasmids_filtered["Plasmid GenBank"].str.split(",")
    genbank_plasmids_filtered_exploded = genbank_plasmids_filtered.explode("Plasmid GenBank")

    genbank_plasmids_filtered_exploded["#Species/genus"] = genbank_plasmids_filtered_exploded[
        "#Species/genus"
    ].str.replace(" ", "_")
    genbank_plasmids_filtered_exploded["#Species/genus"] = genbank_plasmids_filtered_exploded[
        "#Species/genus"
    ].str.replace("/", "_")
    genbank_plasmids_filtered_exploded["#Species/genus"] = genbank_plasmids_filtered_exploded[
        "#Species/genus"
    ].str.replace("'", "")
    genbank_plasmids_filtered_exploded["#Species/genus"] = genbank_plasmids_filtered_exploded[
        "#Species/genus"
    ].str.replace("(", "")
    genbank_plasmids_filtered_exploded["#Species/genus"] = genbank_plasmids_filtered_exploded[
        "#Species/genus"
    ].str.replace(")", "")

    header = ["species", "GBSeq_accession-version"]

    if config["sequences"]:
        custom_fasta_paths = pd.read_csv(
            config["sequences"], sep="\t", header=None, names=["species", "accession", "path"],
        )

        genbank_exploded = genbank_exploded[(~genbank_exploded["species"].isin(custom_fasta_paths["species"]))]
        nuccore_exploded = nuccore_exploded[(~nuccore_exploded["species"].isin(custom_fasta_paths["species"]))]
        assemblies_exploded = assemblies_exploded[(~assemblies_exploded["species"].isin(custom_fasta_paths["species"]))]
        genbank_plasmids_filtered_exploded = genbank_plasmids_filtered_exploded[
            (~genbank_plasmids_filtered_exploded["species"].isin(custom_fasta_paths["species"]))
        ]
        nuccore_plasmids_exploded = nuccore_plasmids_exploded[
            (~nuccore_plasmids_exploded["species"].isin(custom_fasta_paths["species"]))
        ]

    if config["accessions"]:
        custom_accessions = pd.read_csv(config["accessions"], sep="\t", header=None, names=["species", "accession"],)

        genbank_exploded = genbank_exploded[(~genbank_exploded["species"].isin(custom_accessions["species"]))]
        nuccore_exploded = nuccore_exploded[(~nuccore_exploded["species"].isin(custom_accessions["species"]))]
        assemblies_exploded = assemblies_exploded[(~assemblies_exploded["species"].isin(custom_accessions["species"]))]
        genbank_plasmids_filtered_exploded = genbank_plasmids_filtered_exploded[
            (~genbank_plasmids_filtered_exploded["species"].isin(custom_accessions["species"]))
        ]
        nuccore_plasmids_exploded = nuccore_plasmids_exploded[
            (~nuccore_plasmids_exploded["species"].isin(custom_accessions["species"]))
        ]

    genbank_plasmids_to_drop = []
    for key, seq in genbank_plasmids_filtered_exploded.iterrows():
        if key in genbank_plasmids_to_drop:
            continue
        accession = seq["Plasmid GenBank"]
        url = get_accession_ftp_path(accession, config)
        if url != "":
            genbank_plasmids_to_drop.append(key)
    genbank_plasmids_filtered_exploded.drop(genbank_plasmids_to_drop, inplace=True)

    nuccore_plasmids_to_drop = []
    for key, seq in nuccore_plasmids_exploded.iterrows():
        if key in nuccore_plasmids_to_drop:
            continue
        accession = seq["Plasmid RefSeq"]
        url = get_accession_ftp_path(accession, config)
        if url != "":
            nuccore_plasmids_to_drop.append(key)
    nuccore_plasmids_exploded.drop(nuccore_plasmids_to_drop, inplace=True)

    # todo working example the head() needs to leave for a proper run
    genbank_exploded.head(5).to_csv(genbank_genomes_out, sep="\t", header=header, index=False)
    nuccore_exploded.head(5).to_csv(nuccore_genomes_out, sep="\t", header=header, index=False)
    assemblies_exploded.head(5).to_csv(assemblies_out, sep="\t", header=header, index=False)
    genbank_plasmids_filtered_exploded.head(5).to_csv(genbank_plasmids_out, sep="\t", header=header, index=False)
    nuccore_plasmids_exploded.head(5).to_csv(nuccore_plasmids_out, sep="\t", header=header, index=False)


if __name__ == "__main__":
    # redirect all output to the log
    sys.stderr = open(snakemake.log[0], "w")

    entrez_refseq_create_files(
        config=snakemake.config,
        input_file=snakemake.input[0],
        nuccore_genomes_out=snakemake.output[0],
        genbank_genomes_out=snakemake.output[1],
        assemblies_out=snakemake.output[2],
        nuccore_plasmids_out=snakemake.output[3],
        genbank_plasmids_out=snakemake.output[4],
    )
