#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd


def get_total_paths(
    wildcards,
    checkpoints,
    with_entrez_query,
    with_refseq_rep,
    with_custom_sequences,
    with_custom_accessions,
    specific_genera,
):
    """
    Get all the individual fasta file paths for the taxa in our database.
    """

    sequences = pd.DataFrame()

    if with_entrez_query:
        pick_sequences = checkpoints.entrez_pick_sequences.get(query=wildcards.query)
        sequences = pd.read_csv(pick_sequences.output[0], sep="\t")

        if len(sequences) == 0:
            raise RuntimeError("The entrez pick sequences file is empty.")

    if with_refseq_rep:
        refseq_rep_prok = checkpoints.entrez_refseq_accessions.get(
            query=wildcards.query
        )

        refseq_genomes = pd.read_csv(refseq_rep_prok.output.refseq_genomes, sep="\t")
        genbank_genomes = pd.read_csv(refseq_rep_prok.output.genbank_genomes, sep="\t")
        assemblies = pd.read_csv(refseq_rep_prok.output.assemblies, sep="\t")
        refseq_plasmids = pd.read_csv(refseq_rep_prok.output.refseq_plasmids, sep="\t")
        genbank_plasmids = pd.read_csv(
            refseq_rep_prok.output.genbank_plasmids, sep="\t"
        )

        invalid_assemblies = checkpoints.entrez_invalid_assemblies.get(
            query=wildcards.query
        )
        invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep="\t")

        assemblies = assemblies[
            ~assemblies["GBSeq_accession-version"].isin(
                invalid_assembly_sequences["GBSeq_accession-version"]
            )
        ]

        sources = [
            refseq_genomes,
            genbank_genomes,
            assemblies,
            refseq_plasmids,
            genbank_plasmids,
        ]

        if with_entrez_query:
            sources.append(sequences)

        sequences = pd.concat(sources)

    if with_custom_sequences:
        custom_fasta_paths = pd.read_csv(
            config["custom_seq_file"],
            sep="\t",
            header=None,
            names=["species", "GBSeq_accession-version", "path"],
        )

        custom_seqs = custom_fasta_paths[["species", "GBSeq_accession-version"]]

        sequences = sequences.append(custom_seqs)

    if with_custom_accessions:
        custom_accessions = pd.read_csv(
            config["custom_acc_file"],
            sep="\t",
            header=None,
            names=["species", "GBSeq_accession-version"],
        )

        sequences = sequences.append(custom_accessions)

    if specific_genera:
        sequences = sequences[
            sequences["species"].str.contains("|".join(specific_genera))
        ]

    return sequences


def normalise_name(taxon):
    """remove unnecessary characters from a taxon name string."""
    return taxon.replace(" ", "_").replace("[", "").replace("]", "")


def check_unique_taxa_in_custom_input(with_custom_accessions, with_custom_sequences):

    """Checks that custom input files have only one entry per taxon"""

    if with_custom_accessions and with_custom_sequences:
        custom_fasta_paths = pd.read_csv(
            config["custom_seq_file"],
            sep="\t",
            header=None,
            names=["species", "accession", "path"],
        )
        custom_accessions = pd.read_csv(
            config["custom_acc_file"],
            sep="\t",
            header=None,
            names=["species", "accession"],
        )

        taxon_acc = custom_accessions["species"].tolist()
        taxon_seq = custom_fasta_paths["species"].tolist()

        if bool(set(taxon_acc) & set(taxon_seq)):
            raise RuntimeError(
                "You have provided the same taxon both in your custom sequences file and your "
                "custom accessions file. Please pick and keep ONLY one entry from both of these files. "
                "You can only have 1 sequence per chosen taxon in your database."
            )
