#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"

import os
from math import ceil

from psutil import virtual_memory

MEGABYTE = float(1024 ** 2)
MAX_MEM_MB = virtual_memory().total / MEGABYTE

from haystac.workflow.scripts.utilities import get_final_db_paths


def get_db_accessions(_):
    """Get fasta paths in our db"""
    return [
        config["cache"] + f"/ncbi/{orgname}/{accession}.fasta.gz"
        for orgname, accession in get_final_db_paths(checkpoints)
    ]


rule randomise_db_order:
    input:
        get_db_accessions,
    log:
        config["db_output"] + "/bowtie/bt2_random_fasta_paths.log",
    output:
        config["db_output"] + "/bowtie/bt2_random_fasta_paths.txt",
    params:
        seed=config["seed"],
    message:
        "The database genomes are being placed in a random order."
    script:
        "../scripts/random_db_paths.py"


checkpoint calculate_db_chunks:
    input:
        config["db_output"] + "/bowtie/bt2_random_fasta_paths.txt",
    output:
        config["db_output"] + "/bowtie/bt2_idx_chunk_list.txt",
        config["db_output"] + "/bowtie/bt2_idx_chunk_num.txt",
    params:
        query=config["db_output"],
        mem_resources=float(config["mem"]),
        mem_rescale_factor=config["bowtie2_scaling"],
    message:
        "The number of index chunks needed for the filtering alignment are being calculated."
    script:
        "../scripts/calculate_bt2_idx_chunks.py"


rule create_db_chunk:
    input:
        config["db_output"] + "/bowtie/bt2_idx_chunk_list.txt",
    output:
        config["db_output"] + "/bowtie/chunk{chunk_num}.fasta.gz",
    message:
        "Creating chunk {wildcards.chunk_num} of the genome database index."
    threads: 8
    shell:
        "while IFS=$'\t' read -r -a p; "
        "do if [ ${{p[0]}} -eq {wildcards.chunk_num} ]; "
        "then prefix=$(basename ${{p[1]}} .fasta.gz); "
        "bgzip --decompress --stdout --threads 8 ${{p[1]}} | "
        'sed "s/^>/>"${{prefix}}"_/g" | '
        "bgzip --stdout --threads {threads}; fi; done < {input} > {output}"


rule bowtie_index_db_chunk:
    input:
        fasta_chunk=config["db_output"] + "/bowtie/chunk{chunk_num}.fasta.gz",
    log:
        config["db_output"] + "/bowtie/chunk{chunk_num}_index.log",
    output:
        config["db_output"] + "/bowtie/chunk{chunk_num}.1.bt2l",
        config["db_output"] + "/bowtie/chunk{chunk_num}.2.bt2l",
        config["db_output"] + "/bowtie/chunk{chunk_num}.3.bt2l",
        config["db_output"] + "/bowtie/chunk{chunk_num}.4.bt2l",
        config["db_output"] + "/bowtie/chunk{chunk_num}.rev.1.bt2l",
        config["db_output"] + "/bowtie/chunk{chunk_num}.rev.2.bt2l",
    message:
        "Bowtie2 index for chunk {input.fasta_chunk} is being built."
    threads: config["cores"]
    resources:
        mem_mb=lambda wildcards, input: int(
            os.stat(input.fasta_chunk).st_size * config["bowtie2_scaling"] / 1024 ** 2
        ),
        # TODO report memory usage as a function of threads
    conda:
        "../envs/bowtie2.yaml"
    params:
        basename=config["db_output"] + "/bowtie/chunk{chunk_num}",
    priority: 1
    shell:
        "bowtie2-build --large-index --threads {config[cores]} {input.fasta_chunk} {params.basename} &> {log}"


def get_index_db_chunks(_):
    """Get the paths for the index chunks for the filtering bowtie2 alignment"""

    # noinspection PyUnresolvedReferences
    get_chunk_num = checkpoints.calculate_db_chunks.get()
    idx_chunk_total = ceil(float(open(get_chunk_num.output[1]).read().strip()))

    return expand(
        config["db_output"] + "/bowtie/chunk{chunk_num}.1.bt2l",
        chunk_num=[x + 1 if idx_chunk_total > 1 else 1 for x in range(idx_chunk_total)],
    )


rule index_all_db_chunks:
    input:
        get_index_db_chunks,
    output:
        config["db_output"] + "/bowtie/bowtie_index.done",
    message:
        "The bowtie2 indices of the genome database have been built."
    shell:
        "touch {output}"
