#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

##### Target rules #####


rule download_refseq_representative_table:
    output:
        "database_inputs/prok_representative_genomes.txt"
    log:
        "database_inputs/prok_representative_genomes.log"
    benchmark:
        repeat("benchmarks/prok_report_download.benchmark.txt", 3)
    shell:
        "wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_representative_genomes.txt 2> {log}"



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



rule entrez_download_refseq_genbank_sequence:
    output:
        "database/refseq_genbank/{orgname}/{accession}.fasta.gz"
    log:
        "database/refseq_genbank/{orgname}/{accession}.log"
    benchmark:
        repeat("benchmarks/entrez_download_refseq_genbank_sequence_{orgname}_{accession}.benchmark.txt", 3)
    params:
        assembly=False
    script:
        "../scripts/entrez_download_sequence.py"

# ruleorder: entrez_download_refseq_genbank_sequence > entrez_download_sequence



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
        orgname, accession = seq['species'].replace(" ", "_"), seq['GBSeq_accession-version']
        inputs.append('database/refseq_genbank/{orgname}/{accession}.fasta.gz'.format(orgname=orgname,
            accession=accession))

    return inputs



rule entrez_refseq_genbank_multifasta:
    input:
        get_refseq_genome_sequences
    log:
        "{query}/bowtie/{query}_refseq_genbank.log"
    output:
        "{query}/bowtie/{query}_refseq_genbank.fasta.gz"
    benchmark:
        repeat("benchmarks/entrez_refseq_genbank_multifasta_{query}.benchmark.txt", 3)
    shell:
        "cat {input} > {output}"



rule entrez_download_assembly_sequence:
    output:
        "database/refseq_assembly/{orgname}/{accession}.fasta.gz"
    log:
        "database/refseq_assembly/{orgname}/{accession}.log"
    benchmark:
        repeat("benchmarks/entrez_download_assembly_sequence_{orgname}_{accession}.benchmark.txt", 3)
    params:
        assembly=True
    script:
        "../scripts/entrez_download_sequence.py"



# ruleorder: entrez_download_assembly_sequence > entrez_download_sequence



def get_assembly_genome_sequences(wildcards):
    """
    Get all the FASTA sequences for the multi-FASTA file.
    """
    pick_sequences = checkpoints.entrez_refseq_accessions.get(query=wildcards.query)
    assembly_sequences = pd.read_csv(pick_sequences.output[2], sep='\t')

    if len(assembly_sequences) == 0:
        raise RuntimeError("The entrez pick sequences file is empty.")

    inputs = []

    for key, seq in assembly_sequences.iterrows():
        orgname, accession = seq['species'].replace(" ", "_"), seq['GBSeq_accession-version']
        inputs.append('database/refseq_assembly/{orgname}/{accession}.fasta.gz'.format(orgname=orgname,
            accession=accession))

    return inputs



rule entrez_assembly_multifasta:
    input:
        get_assembly_genome_sequences
    log:
        "{query}/bowtie/{query}_assemblies.log"
    output:
        "{query}/bowtie/{query}_assemblies.fasta.gz"
    benchmark:
        repeat("benchmarks/entrez_assembly_multifasta_{query}.benchmark.txt", 3)
    shell:
        "cat {input} > {output}"



rule softlink_refseq_to_database:
    input:
        get_refseq_genome_sequences
    output:
        "database/{query}_softlink_refseq_to_database.done"
    log:
        "database/{query}_softlink_refseq_to_database.log"
    shell:
        "for i in {input}; do ln -sfn `echo $i | cut -d '/' -f 2-3` database; done; touch {output}"



rule softlink_assemblies_to_database:
    input:
        get_assembly_genome_sequences
    output:
        "database/{query}_softlink_assembly_to_database.done"
    log:
        "database/{query}_softlink_assembly_to_database.log"
    shell:
        "for i in {input}; do ln -sfn `echo $i | cut -d '/' -f 2-3` database; done; touch {output}"


rule entrez_refseq_prok_multifasta:
    input:
        assemblies_softlink="database/{query}_softlink_assembly_to_database.done",
        refseq_softlink="database/{query}_softlink_refseq_to_database.done",
        assemblies="{query}/bowtie/{query}_assemblies.fasta.gz",
        refseq="{query}/bowtie/{query}_refseq_genbank.fasta.gz"
    log:
        "{query}/bowtie/{query}_refseq_prok.log"
    output:
        "{query}/bowtie/{query}_refseq_prok.fasta.gz"
    benchmark:
        repeat("benchmarks/entrez_refseq_prok_multifasta_{query}.benchmark.txt", 3)
    shell:
        "cat {input.assemblies} {input.refseq} > {output}"




# todo check the refseq rep files if they're empty

