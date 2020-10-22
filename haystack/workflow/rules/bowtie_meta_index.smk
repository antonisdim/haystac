#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from haystack.workflow.scripts.utilities import normalise_name


rule fasta_idx:
    input:
        config["cache"] + "/ncbi/{orgname}/{accession}.fasta.gz",
    output:
        config["cache"] + "/ncbi/{orgname}/{accession}.fasta.gz.fai",
    benchmark:
        repeat("benchmarks/fasta_idx_{orgname}_{accession}.benchmark.txt", 1)
    message:
        "Indexing fasta file with accession {wildcards.accession} for taxon {wildcards.orgname}."
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"


rule index_database_entrez:
    input:
        config["cache"] + "/ncbi/{orgname}/{accession}.fasta.gz",
    log:
        config["cache"] + "/ncbi/{orgname}/{accession}_index.log",
    output:
        expand(
            config["cache"] + "/ncbi/{orgname}/{accession}.{n}.bt2l", n=[1, 2, 3, 4], allow_missing=True,
        ),
        expand(
            config["cache"] + "/ncbi/{orgname}/{accession}.rev.{n}.bt2l", n=[1, 2], allow_missing=True,
        ),
    benchmark:
        repeat("benchmarks/index_database_{orgname}_{accession}.benchmark.txt", 1)
    message:
        "Preparing the bowtie2 index for genome {wildcards.accession} of taxon {wildcards.orgname}."
    threads: config["bowtie2_threads"]
    params:
        basename=config["cache"] + "/ncbi/{orgname}/{accession}",
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build --large-index --threads {threads} {input} {params.basename} &> {log}"


def get_idx_entrez(_):
    """
    Get all the index paths for the taxa in our database from the entrez query.
    """
    if not config["query"]:
        return []

    # noinspection PyUnresolvedReferences
    pick_sequences = checkpoints.entrez_pick_sequences.get()
    sequences = pd.read_csv(pick_sequences.output[0], sep="\t")

    if len(sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["AccessionVersion"],
        )
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.1.bt2l")
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz.fai")

    return inputs


def get_idx_ref_gen(_):
    """
    Get all the index paths for the taxa in our database from the refseq rep and genbank.
    """
    if not config["refseq_rep"]:
        return []

    # noinspection PyUnresolvedReferences
    refseq_rep_prok = checkpoints.entrez_refseq_accessions.get()

    refseq_genomes = pd.read_csv(refseq_rep_prok.output[0], sep="\t")
    genbank_genomes = pd.read_csv(refseq_rep_prok.output[1], sep="\t")
    refseq_plasmids = pd.read_csv(refseq_rep_prok.output[3], sep="\t")
    genbank_plasmids = pd.read_csv(refseq_rep_prok.output[4], sep="\t")

    sequences = pd.concat([refseq_genomes, genbank_genomes, refseq_plasmids, genbank_plasmids])

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["AccessionVersion"],
        )
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.1.bt2l")
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz.fai")

    return inputs


def get_idx_assembly(_):
    """
    Get all the individual bam file paths for the taxa in our database.
    """
    if not config["refseq_rep"]:
        return []

    # noinspection PyUnresolvedReferences
    refseq_rep_prok = checkpoints.entrez_refseq_accessions.get()

    assemblies = pd.read_csv(refseq_rep_prok.output[2], sep="\t")

    # noinspection PyUnresolvedReferences
    invalid_assemblies = checkpoints.entrez_invalid_assemblies.get()
    invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep="\t")

    assemblies = assemblies[~assemblies["AccessionVersion"].isin(invalid_assembly_sequences["AccessionVersion"])]

    sequences = assemblies

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["AccessionVersion"],
        )
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.1.bt2l")
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz.fai")

    return inputs


def get_idx_custom_seqs():
    """
    Get all the individual bam file paths for the taxa in our database.
    """
    if not config["sequences"]:
        return []

    custom_fasta_paths = pd.read_csv(
        config["sequences"], sep="\t", header=None, names=["species", "accession", "path"],
    )

    sequences = custom_fasta_paths

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["accession"],
        )
        inputs.append(config["cache"] + f"/ncbi/{orgname}/custom_seq-{accession}.1.bt2l")
        inputs.append(config["cache"] + f"/ncbi/{orgname}/custom_seq-{accession}.fasta.gz.fai")

    return inputs


def get_idx_custom_acc():
    """
    Get all the individual bam file paths for the taxa in our database.
    """
    if not config["accessions"]:
        return []

    custom_accessions = pd.read_csv(config["accessions"], sep="\t", header=None, names=["species", "accession"])

    sequences = custom_accessions

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["accession"],
        )
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.1.bt2l")
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz.fai")

    return inputs


rule idx_database:
    input:
        get_idx_entrez,
        get_idx_ref_gen,
        get_idx_assembly,
        get_idx_custom_seqs(),
        get_idx_custom_acc(),
    output:
        config["db_output"] + "/idx_database.done",
    message:
        "All the individual genome indices have been prepared."
    shell:
        "touch {output}"
