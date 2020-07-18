#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evangelos A. Dimopoulos, Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Oxford"
__email__ = "antonisdim41@gmail.com"
__license__ = "MIT"


##### Target rules #####


rule get_sra_fastq_se:
    output:
        temp("sra_data/SE/{accession}.fastq"),
    log:
        temp("sra_data/SE/{accession}.log"),
    params:
        extra="",
    threads: 6
    message:
        "Download SRA file {output} for accession {wildcards.accession}. The log file can be found in {log}."
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump" # TODO we're not using wrappers anywhere else... be consistent!


rule get_sra_fastq_pe:
    output:
        temp("sra_data/PE/{accession}_1.fastq"), # the wildcard name must be accession, pointing to an SRA number
        temp("sra_data/PE/{accession}_2.fastq"),
    log:
        temp("sra_data/PE/{accession}.log"),
    params:
        extra="", # optional extra arguments
    threads: 6 # defaults to 6
    message:
        "Download SRA files {output} for accession {wildcards.accession}. The log file can be found in {log}."
    wrapper:
        "0.51.2/bio/sra-tools/fasterq-dump" # TODO we're not using wrappers anywhere else... be consistent!


rule compress_sra_fastq_se:
    input:
        "sra_data/SE/{accession}.fastq",
    log:
        "sra_data/SE/{accession}_compress.log", # TODO log file is not informative
    output:
        "sra_data/SE/{accession}.fastq.gz",
    message:
        "Compressing the raw fastq file {input} for accession {wildcards.accession} and storing it in {output}. "
        "The log file can be found in {log}."
    shell:
        "gzip -c {input} 1> {output} 2> {log}"


rule compress_sra_fastq_pe:
    input:
        r1="sra_data/PE/{accession}_1.fastq",
        r2="sra_data/PE/{accession}_2.fastq",
    log:
        "sra_data/PE/{accession}_compress.log", # TODO not only is this log file not informative but you immediately overwrite it then you compress r2
    output:
        r1="sra_data/PE/{accession}_R1.fastq.gz",
        r2="sra_data/PE/{accession}_R2.fastq.gz",
    message:
        "Compressing the raw fastq files {input.r1} and {input.r2} for accession {wildcards.accession} "
        "and storing them in {output.r1} and {output.r2}. "
        "The log file can be found in {log}."
    shell:
        "gzip -c {input.r1} 1> {output.r1} 2> {log}; "
        "gzip -c {input.r2} 1> {output.r2} 2> {log}"
