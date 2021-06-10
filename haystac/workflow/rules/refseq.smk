#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

GENOME_REPORTS = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/"
PROK_REFSEQ_REP = "prok_representative_genomes.txt"
EUKARYOTES = "eukaryotes.txt"


rule download_refseq_representative_table:
    output:
        config["db_output"] + "/database_inputs/{table}.txt",
    log:
        config["db_output"] + "/database_inputs/{table}.log",
    message:
        "Downloading the list of representative species from RefSeq."
    conda:
        "../envs/wget.yaml"
    params:
        REFSEQ_REP_URL=GENOME_REPORTS + "{table}.txt",
    shell:
        "wget -O {output} {params.REFSEQ_REP_URL} 2> {log}"


checkpoint entrez_refseq_rep_prok_accessions:
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
    message:
        "Finding if assemblies are not part of the RefSeq database."
    resources:
        entrez_api=1,
    script:
        "../scripts/entrez_invalid_assemblies.py"


checkpoint entrez_refseq_viruses_accessions:
    input:
        config["db_output"] + "/database_inputs/viruses.txt",
    log:
        config["db_output"] + "/entrez/refseq-viruses-seqs.log",
    output:
        refseq_viruses=config["db_output"] + "/entrez/refseq-viruses-seqs.tsv",
    message:
        "Converting the RefSeq table for viruses in smaller table."
    script:
        "../scripts/entrez_refseq_virus_create_files.py"


checkpoint entrez_refseq_eukaryotes_accessions:
    input:
        config["db_output"] + "/database_inputs/eukaryotes.txt",
    log:
        config["db_output"] + "/entrez/refseq-eukaryotes-seqs.log",
    output:
        refseq_euk=config["db_output"] + "/entrez/refseq-eukaryotes-seqs.tsv",
    message:
        "Converting the RefSeq table for eukaryotes in a smaller table."
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/entrez_refseq_eukaryotes_create_files.py"
