#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import pandas as pd

from scripts.rip_utilities import normalise_name

REFSEQ_REP_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_representative_genomes.txt"

##### Target rules #####


rule download_refseq_representative_table:
    output:
        config["db_output"] + "/database_inputs/prok_representative_genomes.txt",
    log:
        config["db_output"] + "/database_inputs/prok_representative_genomes.log",
    benchmark:
        repeat("benchmarks/prok_report_download.benchmark.txt", 3)
    message:
        "Downloading the list of representative species from RefSeq in {output}. "
        "Its log file can be found in {log}."
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
        repeat("benchmarks/entrez_refseq_rep_accessions.benchmark.txt", 3)
    message:
        "Splitting the representative RefSeq table in smaller tables. Table for species whose sequences "
        "can be found in RefSeq:"
        "{output.refseq_genomes}, table for species whose plasmid sequences can be found in RefSeq: "
        "{output.refseq_plasmids}, table for species whose sequences can be found in Genbank: {output.genbank_genomes}, "
        "table for species whose plasmid sequences can be found in Genbank: {output.genbank_plasmids}, "
        "table for species whose sequences can be found in the Assembly database: {output.assemblies}. "
        "The log file can be found in {log}."
    script:
        "../scripts/entrez_refseq_create_files.py"


checkpoint entrez_invalid_assemblies:
    input:
        config["db_output"] + "/entrez/assemblies.tsv",
    log:
        config["db_output"] + "/entrez/invalid-assemblies.log",
    output:
        config["db_output"] + "/entrez/invalid-assemblies.tsv",
    benchmark:
        repeat("benchmarks/entrez_valid_assemblies.benchmark.txt", 1)
    message:
        "Finding assemblies that are not part of the RefSeq database for the accessions in {input}. "
        "The output table can be found in {output} and its log file in {log}."
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_invalid_assemblies.py"


def get_refseq_genome_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_refseq_accessions.get(query=wildcards.query)
    refseq_sequences = pd.read_csv(pick_sequences.output.refseq_genomes, sep="\t")
    genbank_sequences = pd.read_csv(pick_sequences.output.genbank_genomes, sep="\t")
    refseq_plasmids = pd.read_csv(pick_sequences.output.refseq_plasmids, sep="\t")
    genbank_plasmids = pd.read_csv(pick_sequences.output.genbank_plasmids, sep="\t")
    sequences = pd.concat(
        [refseq_sequences, genbank_sequences, refseq_plasmids, genbank_plasmids], axis=0
    )

    if len(sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    inputs = []

    for key, seq in sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["GBSeq_accession-version"],
        )
        inputs.append(
            config["genome_cache_folder"]
            + "/{orgname}/{accession}.fasta.gz".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


rule entrez_refseq_genbank_multifasta:
    input:
        get_refseq_genome_sequences,
    log:
        config["db_output"] + "/bowtie/refseq_genbank.log",
    output:
        config["db_output"] + "/bowtie/refseq_genbank.fasta.gz",
    benchmark:
        repeat("benchmarks/entrez_refseq_genbank_multifasta_.benchmark.txt", 1)
    message:
        "Concatenating all the fasta sequences for all the taxa that can be found in RefSeq and Genbank "
        "in {output}, and its log file can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"


rule entrez_download_assembly_sequence:
    output:
        config["genome_cache_folder"] + "/{orgname}/{accession}.fasta.gz",
    log:
        config["genome_cache_folder"] + "/{orgname}/{accession}.log",
    benchmark:
        repeat(
            "benchmarks/entrez_download_assembly_sequence_{orgname}_{accession}.benchmark.txt",
            1,
        )
    params:
        assembly=True,
    wildcard_constraints:
        accession="[^._]+", # TODO refactor this so we're not reliant on the style of the accession (low priority)
    message:
        "Downloading accession {wildcards.accession} for taxon {wildcards.orgname}. "
        "The downloaded fasta sequence can be found in {output} and its log file in {log}."
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_download_sequence.py"


def get_assembly_genome_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_refseq_accessions.get(query=wildcards.query)
    assembly_sequences = pd.read_csv(pick_sequences.output.assemblies, sep="\t")

    if len(assembly_sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    invalid_assemblies = checkpoints.entrez_invalid_assemblies.get(
        query=wildcards.query
    )
    invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep="\t")

    assembly_sequences = assembly_sequences[
        ~assembly_sequences["GBSeq_accession-version"].isin(
            invalid_assembly_sequences["GBSeq_accession-version"]
        )
    ]

    inputs = []

    for key, seq in assembly_sequences.iterrows():
        orgname, accession = (
            normalise_name(seq["species"]),
            seq["GBSeq_accession-version"],
        )
        inputs.append(
            config["genome_cache_folder"]
            + "/{orgname}/{accession}.fasta.gz".format(
                orgname=orgname, accession=accession
            )
        )

    return inputs


rule entrez_assembly_multifasta:
    input:
        get_assembly_genome_sequences,
    log:
        config["db_output"] + "/bowtie/assemblies.log",
    output:
        config["db_output"] + "/bowtie/assemblies.fasta.gz",
    benchmark:
        repeat("benchmarks/entrez_assembly_multifasta.benchmark.txt", 1)
    message:
        "Concatenating all the fasta sequences for all the taxa that can be found in the Assembly database "
        "in {output}, and its log file can be found in {log}."
    script:
        "../scripts/bowtie2_multifasta.py"


rule entrez_refseq_prok_multifasta:
    input:
        assemblies=config["db_output"] + "/bowtie/assemblies.fasta.gz",
        refseq=config["db_output"] + "/bowtie/refseq_genbank.fasta.gz",
    log:
        config["db_output"] + "/bowtie/refseq_prok.log",
    output:
        config["db_output"] + "/bowtie/refseq_prok.fasta.gz",
    benchmark:
        repeat("benchmarks/entrez_refseq_prok_multifasta.benchmark.txt", 1)
    message:
        "Concatenating input files {input.assemblies} and {input.refseq} in {output}. "
        "Its log file can be found in {log}."
    shell:
        "cat {input.assemblies} {input.refseq} > {output}"
