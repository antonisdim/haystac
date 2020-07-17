#!/usr/bin/env python
# -*- coding: utf-8 -*-

from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider

##### Target rules #####

# TODO why does this file exist? move the contents into `entrez.smk`

ENTREZ_DB_NUCCORE = "nuccore"
ENTREZ_DB_TAXA = "taxonomy"

ENTREZ_RETMODE_XML = "xml"
ENTREZ_RETMODE_TEXT = "text"

ENTREZ_RETTYPE_FASTA = "fasta"
ENTREZ_RETTYPE_GB = "gb"

ENTREZ_RETMAX = 10 ** 9

# noinspection PyUnresolvedReferences
def get_nuccore_xml(wildcards):
    """
    Get all the accession chunks for the {query}-nuccore.tsv file.
    """
    ncbi = NCBIRemoteProvider(
        email=config["entrez"]["email"], keep_local=True, is_default=True
    )
    query = wildcards.query

    accessions = ncbi.search(
        config["entrez"]["queries"][query],
        retmax=ENTREZ_RETMAX,
        retmode=ENTREZ_RETMODE_XML,
        rettype=ENTREZ_RETTYPE_GB,
        idtype="acc",
    )

    input_files = expand("{query}/entrez_alt/{acc}.gb.xml", query=query, acc=accessions)
    input_files.remove("{query}/entrez_alt/NG_052043.1.gb.xml".format(query=query))

    print(input_files)

    # input_files = expand("{acc}.gb.xml", acc=accessions)
    # input_files.remove('NG_052043.1.gb.xml')

    return ncbi.remote(
        input_files, db="nuccore", keep_local=True, force_overwrite=True, batchSize=20
    )


rule entrez_download_xml:
    input:
        get_nuccore_xml,
    log:
        "{query}/entrez_alt/sizes.log",
    output:
        "{query}/entrez_alt/sizes.txt",
    shell:
        "cat {input} > {output}"
