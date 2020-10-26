#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from haystack.workflow.scripts.utilities import normalise_name

REFSEQ_REP_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_representative_genomes.txt"


rule download_refseq_representative_table:
    output:
        config["db_output"] + "/database_inputs/prok_representative_genomes.txt",
    log:
        config["db_output"] + "/database_inputs/prok_representative_genomes.log",
    benchmark:
        repeat("benchmarks/prok_report_download.benchmark.txt", 3)
    message:
        "Downloading the list of representative species from RefSeq."
    conda:
        "../envs/wget.yaml"
    shell:
        "wget -O {output} {REFSEQ_REP_URL} 2> {log}"


checkpoint entrez_refseq_accessions:
    input:
        config["db_output"] + "/database_inputs/prok_representative_genomes.txt",
    log:
        config["db_output"] + "/entrez/refseq-rep-seqs.log",
    output:
        refseq_genomes=config["db_output"] + "/entrez/refseq-genomes.tsv",
        genbank_genomes=config["db_output"] + "/entrez/genbank-genomes.tsv",
        assemblies=config["db_output"] + "/entrez/assemblies.tsv",
        refseq_plasmids=config["db_output"] + "/entrez/refseq-plasmids.tsv",
        genbank_plasmids=config["db_output"] + "/entrez/genbank-plasmids.tsv",
    benchmark:
        repeat("benchmarks/entrez_refseq_rep_accessions.benchmark.txt", 1)
    resources:
        entrez_api=1,
    message:
        "Splitting the representative RefSeq table in smaller tables."
    script:
        "../scripts/entrez_refseq_create_files.py"


checkpoint entrez_invalid_assemblies:
    input:
        config["db_output"] + "/entrez/assemblies.tsv",
    output:
        config["db_output"] + "/entrez/invalid-assemblies.tsv",
    benchmark:
        repeat("benchmarks/entrez_valid_assemblies.benchmark.txt", 1)
    message:
        "Finding if assemblies are not part of the RefSeq database."
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_invalid_assemblies.py"


def get_refseq_genome_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    # noinspection PyUnresolvedReferences
    pick_sequences = checkpoints.entrez_refseq_accessions.get()
    refseq_sequences = pd.read_csv(pick_sequences.output.refseq_genomes, sep="\t")
    genbank_sequences = pd.read_csv(pick_sequences.output.genbank_genomes, sep="\t")
    refseq_plasmids = pd.read_csv(pick_sequences.output.refseq_plasmids, sep="\t")
    genbank_plasmids = pd.read_csv(pick_sequences.output.genbank_plasmids, sep="\t")
    sequences = pd.concat([refseq_sequences, genbank_sequences, refseq_plasmids, genbank_plasmids], axis=0)

    if len(sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["AccessionVersion"],
        )
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz")

    return inputs


rule entrez_refseq_genbank_aggregator:
    input:
        get_refseq_genome_sequences,
    log:
        config["db_output"] + "/bowtie/refseq_genbank.log",
    output:
        config["db_output"] + "/bowtie/refseq_genbank.done",
    benchmark:
        repeat("benchmarks/entrez_refseq_genbank_multifasta_.benchmark.txt", 1)
    message:
        "Aggregating all the fasta sequences for all the taxa that can be found in RefSeq and Genbank."
    shell:
        "touch {output}"


def get_assembly_genome_sequences(_):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    # noinspection PyUnresolvedReferences
    pick_sequences = checkpoints.entrez_refseq_accessions.get()
    assembly_sequences = pd.read_csv(pick_sequences.output.assemblies, sep="\t")

    if len(assembly_sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    if not config["force_accessions"]:
        # noinspection PyUnresolvedReferences
        invalid_assemblies = checkpoints.entrez_invalid_assemblies.get()
        invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep="\t")

        assembly_sequences = assembly_sequences[
            ~assembly_sequences["AccessionVersion"].isin(invalid_assembly_sequences["AccessionVersion"])
        ]

    inputs = []

    for key, seq in assembly_sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["AccessionVersion"],
        )
        inputs.append(config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz")

    return inputs


rule entrez_refseq_assembly_aggregator:
    input:
        get_assembly_genome_sequences,
    log:
        config["db_output"] + "/bowtie/assemblies.log",
    output:
        config["db_output"] + "/bowtie/assemblies.done",
    benchmark:
        repeat("benchmarks/entrez_assembly_multifasta.benchmark.txt", 1)
    message:
        "Aggregating all the fasta sequences for all the taxa that can be found in the assembly database."
    shell:
        "touch {output}"


rule entrez_refseq_prok_aggregator:
    input:
        assemblies=config["db_output"] + "/bowtie/assemblies.done",
        refseq=config["db_output"] + "/bowtie/refseq_genbank.done",
    log:
        config["db_output"] + "/bowtie/refseq_prok.log",
    output:
        config["db_output"] + "/bowtie/refseq_prok.done",
    benchmark:
        repeat("benchmarks/entrez_refseq_prok_multifasta.benchmark.txt", 1)
    message:
        "Aggregating input files {input.assemblies} and {input.refseq}."
    shell:
        "touch {output}"
