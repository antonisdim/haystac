#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

from itertools import chain

from haystack.workflow.scripts.utilities import get_total_paths


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


rule idx_database:
    input:
        list(
            chain.from_iterable(
                (
                    config["cache"] + f"/ncbi/{orgname}/{accession}.1.bt2l",
                    config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz.fai",
                )
                for orgname, accession in get_total_paths(checkpoints, config)
            )
        ),
    output:
        config["db_output"] + "/idx_database.done",
    message:
        "All the individual genome indices have been prepared."
    shell:
        "touch {output}"
