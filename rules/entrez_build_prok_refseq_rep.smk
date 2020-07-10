#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

REFSEQ_REP_URL = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_representative_genomes.txt'

##### Target rules #####


rule download_refseq_representative_table:
    output:
        "database_inputs/prok_representative_genomes.txt"
    log:
        "database_inputs/prok_representative_genomes.log"
    benchmark:
        repeat("benchmarks/prok_report_download.benchmark.txt", 3)
    shell:
        "wget -O {output} {REFSEQ_REP_URL} 2> {log}"



checkpoint entrez_refseq_accessions:
    input:
        "database_inputs/prok_representative_genomes.txt"
    log:
        "{query}/entrez/{query}-refseq-seqs.log"
    output:
        refseq_genomes="{query}/entrez/{query}-refseq-genomes.tsv",
        genbank_genomes="{query}/entrez/{query}-genbank-genomes.tsv",
        assemblies="{query}/entrez/{query}-assemblies.tsv",
        refseq_plasmids="{query}/entrez/{query}-refseq-plasmids.tsv",
        genbank_plasmids="{query}/entrez/{query}-genbank-plasmids.tsv"
    benchmark:
        repeat("benchmarks/entrez_refseq_accessions_{query}.benchmark.txt", 3)
    script:
        "../scripts/entrez_refseq_create_files.py"



checkpoint entrez_invalid_assemblies:
    input:
        "{query}/entrez/{query}-assemblies.tsv"
    log:
        "{query}/entrez/{query}-invalid-assemblies.log"
    output:
        "{query}/entrez/{query}-invalid-assemblies.tsv"
    benchmark:
        repeat("benchmarks/entrez_valid_assemblies_{query}.benchmark.txt", 1)
    script:
        "../scripts/entrez_invalid_assemblies.py"



def get_refseq_genome_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_refseq_accessions.get(query=wildcards.query)
    refseq_sequences = pd.read_csv(pick_sequences.output[0], sep='\t')
    genbank_sequences = pd.read_csv(pick_sequences.output[1], sep='\t')
    refseq_plasmids = pd.read_csv(pick_sequences.output[3], sep='\t')
    genbank_plasmids = pd.read_csv(pick_sequences.output[4], sep='\t')
    sequences = pd.concat([refseq_sequences, genbank_sequences, refseq_plasmids, genbank_plasmids], axis=0)

    if len(sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    inputs = []

    for key, seq in sequences.iterrows():
        # TODO move the .replace() code into a function called normalize_name() and update any other instances
        orgname, accession = seq['species'].replace(" ", "_").replace("[", "").replace("]", ""), seq['GBSeq_accession-version']
        inputs.append('database/{orgname}/{accession}.fasta.gz'.format(orgname=orgname, accession=accession))

    return inputs



rule entrez_refseq_genbank_multifasta:
    input:
        get_refseq_genome_sequences
    log:
        "{query}/bowtie/{query}_refseq_genbank.log"
    output:
        "{query}/bowtie/{query}_refseq_genbank.fasta.gz"
    benchmark:
        repeat("benchmarks/entrez_refseq_genbank_multifasta_{query}.benchmark.txt", 1)
    script:
        "../scripts/bowtie2_multifasta.py"



rule entrez_download_assembly_sequence:
    output:
        "database/{orgname}/{accession}.fasta.gz"
    log:
        "database/{orgname}/{accession}.log"
    benchmark:
        repeat("benchmarks/entrez_download_assembly_sequence_{orgname}_{accession}.benchmark.txt", 1)
    params:
        assembly=True
    wildcard_constraints:
        # TODO refactor this so we're not reliant on the style of the accession (low priority)
        accession="[^._]+"
    resources:
        # TODO add this to every other rule that needs it
        entrez_api=1
    script:
        "../scripts/entrez_download_sequence.py"



def get_assembly_genome_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_refseq_accessions.get(query=wildcards.query)
    assembly_sequences = pd.read_csv(pick_sequences.output[2], sep='\t')

    if len(assembly_sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    invalid_assemblies = checkpoints.entrez_invalid_assemblies.get(query=wildcards.query)
    invalid_assembly_sequences = pd.read_csv(invalid_assemblies.output[0], sep='\t')

    assembly_sequences = assembly_sequences[~assembly_sequences['GBSeq_accession-version'].isin(invalid_assembly_sequences['GBSeq_accession-version'])]

    inputs = []

    for key, seq in assembly_sequences.iterrows():
        orgname, accession = seq['species'].replace(" ", "_").replace("[", "").replace("]", ""), seq['GBSeq_accession-version']
        inputs.append('database/{orgname}/{accession}.fasta.gz'.format(orgname=orgname, accession=accession))

    return inputs



rule entrez_assembly_multifasta:
    input:
        get_assembly_genome_sequences
    log:
        "{query}/bowtie/{query}_assemblies.log"
    output:
        "{query}/bowtie/{query}_assemblies.fasta.gz"
    benchmark:
        repeat("benchmarks/entrez_assembly_multifasta_{query}.benchmark.txt", 1)
    script:
        "../scripts/bowtie2_multifasta.py"



rule entrez_refseq_prok_multifasta:
    input:
        assemblies="{query}/bowtie/{query}_assemblies.fasta.gz",
        refseq="{query}/bowtie/{query}_refseq_genbank.fasta.gz"
    log:
        "{query}/bowtie/{query}_refseq_prok.log"
    output:
        "{query}/bowtie/{query}_refseq_prok.fasta.gz"
    benchmark:
        repeat("benchmarks/entrez_refseq_prok_multifasta_{query}.benchmark.txt", 1)
    shell:
        "cat {input.assemblies} {input.refseq} > {output}"



