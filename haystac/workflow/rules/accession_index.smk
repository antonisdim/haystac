#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
from itertools import chain

from haystac.workflow.scripts.utilities import get_total_paths


rule samtools_index_accession:
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


rule bowtie_index_accession:
    input:
        config["cache"] + "/ncbi/{orgname}/{accession}.fasta.gz",
    log:
        config["cache"] + "/ncbi/{orgname}/{accession}_index.log",
    output:
        config["cache"] + "/ncbi/{orgname}/{accession}.1.bt2l",
        config["cache"] + "/ncbi/{orgname}/{accession}.2.bt2l",
        config["cache"] + "/ncbi/{orgname}/{accession}.3.bt2l",
        config["cache"] + "/ncbi/{orgname}/{accession}.4.bt2l",
        config["cache"] + "/ncbi/{orgname}/{accession}.rev.1.bt2l",
        config["cache"] + "/ncbi/{orgname}/{accession}.rev.2.bt2l",
    benchmark:
        repeat("benchmarks/index_database_{orgname}_{accession}.benchmark.txt", 1)
    message:
        "Preparing the bowtie2 index for genome {wildcards.accession} of taxon {wildcards.orgname}."
    threads: lambda wildcards, input: 1 if os.stat(input[0]).st_size / (1024 ** 2) < 100 else config["cores"]
    resources:
        mem_mb=lambda wildcards, input: int(os.stat(input[0]).st_size * config["bowtie2_scaling"] / 1024 ** 2), # TODO report memory usage as a function of threads
    params:
        basename=config["cache"] + "/ncbi/{orgname}/{accession}",
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build --large-index --threads {threads} {input} {params.basename} &> {log}"


def get_index_paths(_):
    """Get paths for db indices"""
    return list(
        chain.from_iterable(
            (
                config["cache"] + f"/ncbi/{orgname}/{accession}.1.bt2l",
                config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz.fai",
            )
            for orgname, accession in get_total_paths(checkpoints, config)
        )
    )


rule index_all_accessions:
    input:
        get_index_paths,
    output:
        config["db_output"] + "/idx_database.done",
    message:
        "All the individual genome indices have been prepared."
    shell:
        "touch {output}"
