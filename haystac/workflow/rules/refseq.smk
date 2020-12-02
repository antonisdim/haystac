#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

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
